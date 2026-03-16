#include "fem_matrix.h"

#include "P1.h"
// 解释FEMatrix类的两个虚函数:mvp和sum 他们根据fem_tpye的不同,调用对应函数
void FEMatrix::mvp(const double *x, double *y) const //矩阵向量乘法y=Ax,y是输出
{
	switch (fem_type) { 
	case FEMatrix::P1_cst:
		mvp_P1_cst(*this, x, y);
		return;
	case FEMatrix::P1_sym:
		mvp_P1_sym(*this, x, y);
		return;
	case FEMatrix::P1_gen:
		mvp_P1_gen(*this, x, y);
		return;
	}
}

double FEMatrix::sum() const //计算矩阵元素和
{
	switch (fem_type) {
	case FEMatrix::P1_cst:
		return sum_P1_cst(*this);
	case FEMatrix::P1_sym:
		return sum_P1_sym(*this);
	case FEMatrix::P1_gen:
		return sum_P1_gen(*this);
	default:
		return 0;
	}
}

/* switch 是多分支选择语句
 * 当嵌套的if 比较少时（三个以内），用 if 编写程序会比较简洁。但是当选择的分支比较多时，用 switch 语句来处理多分支选择。
 * 		switch(表达式)
		{
		case 常量表达式1:语句1;
		case 常量表达式2:语句2;
		...
		case 常量表达式n:语句n;
		default:语句;
		}            
 *
 * switch() 括号内的“表达式”必须是整数类型或者可以转换为整型的数值类型。比如：byte、short、int、char、也可以直接是整数或字符常量还有枚举，哪怕是负数都可以。
 * 但是float、double、long和String类型是不能作用在switch语句上的。
 * 执行完一个case后面的语句后，流程控制转移到下一个case继续执行。如果你只想执行这一个case语句，不想执行其他case，那么就需要在这个case语句后面加上break，跳出switch语句。
 */