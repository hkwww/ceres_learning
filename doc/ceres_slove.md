# ceres中slove

### slove过程中重要的函数

- **Slover::solve**

```c++
void Solver::Solve(const Solver::Options& options,
​                   Problem* problem,
​                   Solver::Summary* summary)

```

```
void Solver::Solve(const Solver::Options& options,
​                   Problem* problem,
​                   Solver::Summary* summary) {
  using internal::PreprocessedProblem;
  using internal::Preprocessor;
  using internal::ProblemImpl;
  using internal::Program;
  using internal::scoped_ptr;
  using internal::WallTimeInSeconds;

  CHECK_NOTNULL(problem);
  CHECK_NOTNULL(summary);

  double start_time = WallTimeInSeconds();
  *summary = Summary();
  if (!options.IsValid(&summary->message)) {
    LOG(ERROR) << "Terminating: " << summary->message;
    return;
  }

  ProblemImpl* problem_impl = problem->problem_impl_.get();
  Program* program = problem_impl->mutable_program();
  PreSolveSummarize(options, problem_impl, summary);

  // 设置线程池的线程数
  problem_impl->context()->EnsureMinimumThreads(options.num_threads - 1);

  program->SetParameterBlockStatePtrsToUserStatePtrs();

  //设置梯度检测
  scoped_ptr<internal::ProblemImpl> gradient_checking_problem;
  internal::GradientCheckingIterationCallback gradient_checking_callback;
  Solver::Options modified_options = options;
  if (options.check_gradients) {
    modified_options.callbacks.push_back(&gradient_checking_callback);
    gradient_checking_problem.reset(
        CreateGradientCheckingProblemImpl(
            problem_impl,
            options.gradient_check_numeric_derivative_relative_step_size,
            options.gradient_check_relative_precision,
            &gradient_checking_callback));
    problem_impl = gradient_checking_problem.get();
    program = problem_impl->mutable_program();
  }

 //预处理,根据算法分为(LineSearchPreprocessor::Preprocess() && TrustRegionPreprocessor::Preprocess())
 //初始化线性求解器;初始化evaluator;初始化InnerIterationMinimizer;初始化最小化option
  scoped_ptr<Preprocessor> preprocessor(
      Preprocessor::Create(modified_options.minimizer_type));
  PreprocessedProblem pp;
  const bool status = preprocessor->Preprocess(modified_options, problem_impl, &pp);

  //有关舒尔补的操作
  if (IsSchurType(pp.linear_solver_options.type)) {
    int row_block_size;
    int e_block_size;
    int f_block_size;
    DetectStructure(*static_cast<internal::BlockSparseMatrix*>(
                    pp.minimizer_options.jacobian.get())
                    ->block_structure(),
                    pp.linear_solver_options.elimination_groups[0],
                    &row_block_size,
                    &e_block_size,
                    &f_block_size);
    summary->schur_structure_given =
        SchurStructureToString(row_block_size, e_block_size, f_block_size);
    internal::GetBestSchurTemplateSpecialization(&row_block_size,
                                                 &e_block_size,
                                                 &f_block_size);
    summary->schur_structure_used =
        SchurStructureToString(row_block_size, e_block_size, f_block_size);
  }


  summary->fixed_cost = pp.fixed_cost; 
  //预处理的时间
  summary->preprocessor_time_in_seconds = WallTimeInSeconds() - start_time;

  //循环迭代主要在 Minimize(&pp, summary)函数中
  if (status) {
    const double minimizer_start_time = WallTimeInSeconds();
    Minimize(&pp, summary);
    summary->minimizer_time_in_seconds =
       WallTimeInSeconds() - minimizer_start_time;
  } else {
    summary->message = pp.error;
  }

  //后处理部分
  const double postprocessor_start_time = WallTimeInSeconds();
  problem_impl = problem->problem_impl_.get();
  program = problem_impl->mutable_program();

  program->SetParameterBlockStatePtrsToUserStatePtrs();
  program->SetParameterOffsetsAndIndex();
  PostSolveSummarize(pp, summary);
  summary->postprocessor_time_in_seconds =
      WallTimeInSeconds() - postprocessor_start_time;


  if (gradient_checking_callback.gradient_error_detected()) {
    summary->termination_type = FAILURE;
    summary->message = gradient_checking_callback.error_log();
  }

  summary->total_time_in_seconds = WallTimeInSeconds() - start_time;

}
```







- **Minimize** 

  解决最小二乘问题
  
  ```
  void Minimize(internal::PreprocessedProblem* pp,
                Solver::Summary* summary) {
    using internal::Program;
    using internal::scoped_ptr;
    using internal::Minimizer;
  
    Program* program = pp->reduced_program.get();
  
    if (pp->reduced_program->NumParameterBlocks() == 0) {
      summary->message = "Function tolerance reached. "
          "No non-constant parameter blocks found.";
      summary->termination_type = CONVERGENCE;
      VLOG_IF(1, pp->options.logging_type != SILENT) << summary->message;
      summary->initial_cost = summary->fixed_cost;
      summary->final_cost = summary->fixed_cost;
      return;
    }
  
    const Vector original_reduced_parameters = pp->reduced_parameters;
    //解最小二乘问题,根据优化策略派生为(LineSearchMinimizer::Minimize() &&  TrustRegionMinimizer::Minimize())
    scoped_ptr<Minimizer> minimizer(
        Minimizer::Create(pp->options.minimizer_type));
    minimizer->Minimize(pp->minimizer_options,
                        pp->reduced_parameters.data(),
                        summary);
  
    program->StateVectorToParameterBlocks(
        summary->IsSolutionUsable()
        ? pp->reduced_parameters.data()
        : original_reduced_parameters.data());
    program->CopyParameterBlockStateToUserState();
  
  }
  ```



    void TrustRegionMinimizer::Minimize(const Minimizer::Options& options,
                                        double* parameters,
                                        Solver::Summary* solver_summary) {
      start_time_in_secs_ = WallTimeInSeconds();
      iteration_start_time_in_secs_ = start_time_in_secs_;
      Init(options, parameters, solver_summary);
      RETURN_IF_ERROR_AND_LOG(IterationZero());
    
      step_evaluator_.reset(new TrustRegionStepEvaluator(
          x_cost_,
          options_.use_nonmonotonic_steps
              ? options_.max_consecutive_nonmonotonic_steps
              : 0));
    
      //该函数决定是否停止迭代,具体看下一个代码块
      while (FinalizeIterationAndCheckIfMinimizerCanContinue()) {
        //重新构建interaction_summary
        iteration_start_time_in_secs_ = WallTimeInSeconds();
        iteration_summary_ = IterationSummary();
        iteration_summary_.iteration =
            solver_summary->iterations.back().iteration + 1;
            
        //主要运算函数,计算trust regin还有雅克比等
        RETURN_IF_ERROR_AND_LOG(ComputeTrustRegionStep());
        if (!iteration_summary_.step_is_valid) {
          RETURN_IF_ERROR_AND_LOG(HandleInvalidStep());
          continue;
        }
    
        //最小二乘问题是否有上下界
        if (options_.is_constrained) {
          // Use a projected line search to enforce the bounds constraints
          // and improve the quality of the step.
          DoLineSearch(x_, gradient_, x_cost_, &delta_);
        }
    
        ComputeCandidatePointAndEvaluateCost();
        DoInnerIterationsIfNeeded();
    
        if (ParameterToleranceReached()) {
          return;
        }
    
        if (FunctionToleranceReached()) {
          return;
        }
        if (IsStepSuccessful()) {
          RETURN_IF_ERROR_AND_LOG(HandleSuccessfulStep());
          continue;
        }
    
        HandleUnsuccessfulStep();
        }
        }



- **TrustRegionMinimizer::FinalizeIterationAndCheckIfMinimizerCanContinue**

    该函数用来判断是否停止迭代,主要有三个作用:

    1. 给interaction summary添加时间信息.

    2. callback的调用.

    3. 检测迭代是否停止:

       a. 工作时间

       b. 迭代次数

       c. 梯度的模长

       d. 信赖域的大小



      bool TrustRegionMinimizer::FinalizeIterationAndCheckIfMinimizerCanContinue() {
      //对interaction summary进行赋值
        if (iteration_summary_.step_is_successful) {
          ++solver_summary_->num_successful_steps;
          if (x_cost_ < minimum_cost_) {
            minimum_cost_ = x_cost_;
            VectorRef(parameters_, num_parameters_) = x_;
            iteration_summary_.step_is_nonmonotonic = false;
          } else {
            iteration_summary_.step_is_nonmonotonic = true;
          }
        } else {
          ++solver_summary_->num_unsuccessful_steps;
        }
        iteration_summary_.trust_region_radius = strategy_->Radius();
        iteration_summary_.iteration_time_in_seconds =
            WallTimeInSeconds() - iteration_start_time_in_secs_;
        iteration_summary_.cumulative_time_in_seconds =
            WallTimeInSeconds() - start_time_in_secs_ +
            solver_summary_->preprocessor_time_in_seconds;
    
        solver_summary_->iterations.push_back(iteration_summary_);
    
        //处理callback函数
        if (!RunCallbacks(options_, iteration_summary_, solver_summary_)) {
          return false;
        }
    
        //检测求解时间是否超过设定的最长时间
        if (MaxSolverTimeReached()) {
          return false;
        }
    
        //求解迭代次数是否超过最大迭代次数
        if (MaxSolverIterationsReached()) {
          return false;
        }
    
        //梯度的最大范数是否小于设定值
        if (GradientToleranceReached()) {
          return false;
        }
    
        //trust regin的大小是否小于设定值
        if (MinTrustRegionRadiusReached()) {
          return false;
         }
        return true;
        }

  

