



# Ceres库中构建problem的重要类

## Class Problem

Problem类应该是ceres中的核心.

### Members Struct

**Options**

​	结构体,优化的选项.

**EvaluateOptions**

​	结构体,Evaluate函数的选项.



### Members

**internal::scoped_ptr< internal::ProblemImpl > problem_impl_**

​	Class problem只是个接口,主要功能函数都在problem_impl_(Class ProblemImpl)完成.

### Important Member Function

**AddResidualBlock**

​	添加残差块添加到problem中,有多种形式,主要为:

```c++
  ResidualBlockId AddResidualBlock(
      CostFunction* cost_function,
      LossFunction* loss_function,
      const std::vector<double*>& parameter_blocks);
```

代价函数:包含了参数模块的维度信息，内部使用仿函数定义误差函数的计算方式.

损失函数:处理参数中有野值的情况.

参数模块:待优化的参数(隐式传递参数)

**AddParameterBlock**

​	添加参数块(重复调用相同参数会被无视)

```c++
 void AddParameterBlock(double* values,
                        int size,
                        LocalParameterization* local_parameterization);
```

LocalParameterization类的作用是解决非线性优化中的过参数化问题。所谓过参数化，即待优化参数的实际自由度小于参数本身的自由度。例如在SLAM中，当采用四元数表示位姿时，由于四元数本身的约束（模长为1），实际的自由度为3而非4。此时，若直接传递四元数进行优化，冗余的维数会带来计算资源的浪费，需要使用Ceres预先定义的QuaternionParameterization对优化参数进行重构.
LocalParameterization见下面.

**SetParameterBlockConstant**

​	设置参数块为不变的,优化过程中不改变该变量,Class ParameterBlock中有一个  is_constant_ 参数,设置为true.

**SetParameterBlockVariable**

​	设置参数块为可变的,优化过程中优化的量,Class ParameterBlock中有一个  is_constant_ 参数,设置为false.

**Evaluate**

​	该函数紧跟在参数赋值后，在给定的参数位置求解Problem，给出当前位置处的cost、梯度以及Jacobian矩阵.





## class LocalParameterization

该类为虚基类,ceres定义几种类集成了该类,也可以自定义.

```c++
class CERES_EXPORT LocalParameterization {
 public:
  virtual ~LocalParameterization();
  virtual bool Plus(const double* x,
                    const double* delta,
                    double* x_plus_delta) const = 0;

  virtual bool ComputeJacobian(const double* x, double* jacobian) const = 0;

  virtual bool MultiplyByJacobian(const double* x,
                                  const int num_rows,
                                  const double* global_matrix,
                                  double* local_matrix) const;

  virtual int GlobalSize() const = 0;

  virtual int LocalSize() const = 0;
};


```

若是自定义的函数,上述成员函数中,需要我们改写的是:

**GlobalSize()** 本身的参数的个数(例如QuaternionParameterization ,GlobalSize()  = 4,四元数本身上是四维的变量)

 **LocalSize()** 实际优化的参数个数(例如QuaternionParameterization ,LocalSize()  = 3,ceres采用旋转矢量就行优化,维数为3)

**Plus() ComputeJacobian()** 这需要根据不同参数来定,比如q元素相加





## Class CostFunction

CostFunction用来计算残差和雅克比矩阵,ceres提供不同类型的CostFunction,都是从这个类开始继承的.

```c++
class CostFunction {
 public:
  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) = 0;
  const vector<int32>& parameter_block_sizes();
  int num_residuals() const;

 protected:
  vector<int32>* mutable_parameter_block_sizes();
  void set_num_residuals(int num_residuals);
 private:
  std::vector<int32> parameter_block_sizes_;
  int num_residuals_;
};
```

Evaluate()函数要根据参数计算残差和雅克比.



## Class SizedCostFunction

Class SizedCostFunction继承了Class CostFunction

当在编译前已知优化参数的大小,残差的大小的时候可以使用这个类,只需要重写一下Evaluate函数就好了

```c++
template<int kNumResiduals,
         int N0 = 0, int N1 = 0, int N2 = 0, int N3 = 0, int N4 = 0,
         int N5 = 0, int N6 = 0, int N7 = 0, int N8 = 0, int N9 = 0>
class SizedCostFunction : public CostFunction {
 public:
  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const = 0;
};
```





## class AutoDiffCostFunction

自动求导的代价函数,继承了Class SizedCostFunction.

```c++
template <typename CostFunctor,
       int kNumResiduals,  // Number of residuals, or ceres::DYNAMIC.
       int N0,       // Number of parameters in block 0.
       int N1 = 0,   // Number of parameters in block 1.
       int N2 = 0,   // Number of parameters in block 2.
       int N3 = 0,   // Number of parameters in block 3.
       int N4 = 0,   // Number of parameters in block 4.
       int N5 = 0,   // Number of parameters in block 5.
       int N6 = 0,   // Number of parameters in block 6.
       int N7 = 0,   // Number of parameters in block 7.
       int N8 = 0,   // Number of parameters in block 8.
       int N9 = 0>   // Number of parameters in block 9.
class AutoDiffCostFunction : public
SizedCostFunction<kNumResiduals, N0, N1, N2, N3, N4, N5, N6, N7, N8, N9> {
 public:
  explicit AutoDiffCostFunction(CostFunctor* functor);
  // Ignore the template parameter kNumResiduals and use
  // num_residuals instead.
  AutoDiffCostFunction(CostFunctor* functor, int num_residuals);
};
```

需要自定义一个costfunction类作为类模板的CostFunctor,重写一下operator()函数,要用类模板.

使用说明:

```c++
CostFunction* cost_function
    = new AutoDiffCostFunction<MyScalarCostFunctor, 1, 2, 2>(
        new MyScalarCostFunctor(1.0));              ^  ^  ^
                                                    |  |  |
                        Dimension of residual ------+  |  |
                        Dimension of x ----------------+  |
                        Dimension of y -------------------+
```





## Class DynamicAutoDiffCostFunction

Class AutoDiffCostFunction适用于事先知道参数的大小和数量,且不超过10个参数块,但在许多问题上不适用,比如曲线拟合,所以需要Class DynamicAutoDiffCostFunction.

Class DynamicAutoDiffCostFunction继承了虚基类Class DynamicCostFunction,该类是从Class CostFunction继承过来的.

使用时需要定义一个类作为类模板CostFunctor,重写一下operator()函数

```c++
struct MyCostFunctor {
  template<typename T>
  bool operator()(T const* const* parameters, T* residuals) const {
  }
}
```

之后还要调用AddParameterBlock()和SetNumResiduals()函数来动态添加参数块.





## Class NumericDiffCostFunction

在有些时候无法编写出数值上的CostFunctor,比如需要调用库函数的时候,此时就需要Class NumericDiffCostFunction.

像auto diff是通过公式自动求导,而NumericDiff则是通过数值差来求导,有三种求导方式:

**FORWARD**  使用该公式来求导,速度最快但经度最低
$$
f′(x)=\frac {f(x+h)−f(x)}h
$$
**CENTRAL** 默认使用的方法
$$
f′(x)=\frac {f(x+h)−f(x-h)}{2h}
$$
**RIDDERS** 这种方法是在不同尺度执行中心差来估计导数,是种自适应的方法,经度很高,但是会造成计算的成本.