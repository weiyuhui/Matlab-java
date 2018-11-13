package test;

//import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;

import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NonMonotonicSequenceException;
import org.apache.commons.math3.exception.NotFiniteNumberException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;

import flanagan.interpolation.CubicInterpolation;
import flanagan.interpolation.CubicSpline;
import flanagan.interpolation.PolyCubicSpline;
import flanagan.interpolation.PolylineSimplification;
import flanagan.io.Db;
import flanagan.math.Point;
import interpolation.Interpolation;
import interpolation.Spline;

public class TestApcheInterspline {

	public static void main(String[] args) {

		/*
		 * 方法1：Apache的工具包
		 * 小数点后第1或2位的差别
		 * 
		 */

//		double[] x = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
//		double[] y = { 0.0, 0.8414, 0.9092, 0.1411, -0.7568, -0.9589, -0.2794, 0.6569, 0.9893, 0.4121, -0.5440 };
//		double[] xx = { 0, 0.2500, 0.500, 0.7500, 1, 1.2500, 1.500, 1.7500, 2, 2.2500, };
//
//		// 样条计算法
//		SplineInterpolator sp = new SplineInterpolator();
//		PolynomialSplineFunction f = sp.interpolate(x, y);
//		for (int i = 0; i < xx.length; i++) {
//			System.out.println(f.value(xx[i]));
//		}
//		
//		float[] x2 = {(float)0.0, (float)1.0, (float)2.0, (float)3.0, (float)4.0};
//		float[] y2 = { (float)0.0, (float)0.8414, (float)0.9092, (float)0.1411, (float)-0.7568};
//		float[] xx2 = { 0, (float)0.2500, (float)0.500, (float)0.7500, 1, (float)1.2500, (float)1.500, 
//				(float)1.7500, 2, (float)2.2500, };
//		
//		Spline f2 = Spline.createMonotoneCubicSpline(x2, y2);//单调的X和Y
//		for (int i = 0; i < xx2.length; i++) {
//			System.out.println(f2.interpolate(xx2[i]));
//		}

		/*
		 * 方法2：
		 * Michael T Flanagan的工具包 CubicSpline 小数点后第三位的差别。
		 * 
		 */
		double[] wavelength = { 1,8,15,23,31,39,47,55,63,71,79,87,94,102,110,118,126,134,42,149,157,165,173,181,188,196,204,212,220,227,235 };
		double[] refrindex = {1.34026,1.3402,1.2658,1.2500,1.2822,1.2658,1.2658,1.2500,1.2822,1.2658,1.2658,1.2658,1.2658,1.28205,1.2658,1.2658,1.2987,1.2658,1.2658,1.2658,1.2822,1.2822,1.2987,1.2822,1.2987,1.2500,1.2821,1.2821,1.2821,1.2821,1.2658};
		
		double x1, y1;
		x1 = 5.7;
		
		/*
		 * 比较接近的方法了，小数点后第三位的差别
		 * 结果如数据本身有关系
		 * 是否与数据长度有关系呢？
		 * 测试更长的数据：
		 */
		int order = 3;
		double y = poly_interpolate(wavelength, refrindex, 1, 3);
		System.out.println(y);
		System.out.println(poly_interpolate(wavelength, refrindex, 2, order));
		System.out.println(poly_interpolate(wavelength, refrindex, 3, order));
		System.out.println(poly_interpolate(wavelength, refrindex, 4, order));
		System.out.println(poly_interpolate(wavelength, refrindex, 5, order));
		System.out.println(poly_interpolate(wavelength, refrindex, 6, order));
		System.out.println(poly_interpolate(wavelength, refrindex, 7, order));
		System.out.println(poly_interpolate(wavelength, refrindex, 8, order));
		System.out.println(poly_interpolate(wavelength, refrindex, 9, order));
		System.out.println(poly_interpolate(wavelength, refrindex, 10, order));
		


//		System.out.println("The refractive index of fused quartz at " + x1 * 1.0e9 + " nm is " + y1);
//		x1 = 5.9e-7;
//		y1 = cs.interpolate(x1);
//		System.out.println("The refractive index of fused quartz at " + x1 * 1.0e9 + " nm is " + y1);

		// CubicSpline cs =new CubicSpline(x, y);
		// for(int i = 0; i<xx.length; i++){
		// System.out.println(cs.interpolate(xx[i]));
		// }
		
		/*
		 * Q1:约束系数如何确定的
		 */
		
//		 // Array of wavelengths (m)
//        double[] wavelength = {185.0e-9, 214.0e-9, 275.0e-9, 361.0e-9, 509.0e-9, 589.0e-9, 656.0e-9};
//        // Array of corresponding refractive indices
//        double[] refrindex = {1.57464, 1.53386, 1.49634, 1.47503, 1.4619, 1.4583, 1.4564};
//
//        // Create a CubicInterpolation instance and initialise it to the data stored in the arrays wavelength and refrindex
//        // Numerical differencing on entered data
//        CubicInterpolation ci1 =new CubicInterpolation(wavelength, refrindex, 0);
//
//        // Create a CubicInterpolation instance and initialise it to the data stored in the arrays wavelength and refrindex
//        // Numerical differencing on interpolated data
//        CubicInterpolation ci2 =new CubicInterpolation(wavelength, refrindex, 1);
//
//        // Create a CubicSpline instance and initialise it to the data stored in the arrays wavelength and refrindex
//        CubicSpline cs =new CubicSpline(wavelength, refrindex);
//
//        // First interpolation at a wavelength of 250 nm
//        double x1 = 2.5e-7;
//
//        double y1=ci1.interpolate(x1);
//        System.out.println("Cubic interpolation (numerical differencing on entered data):");
//        System.out.println("The refractive index of fused quartz at " + x1*1.0e9 + " nm is "+ y1);
//
//        double y2=ci2.interpolate(x1);
//        System.out.println("Cubic interpolation (numerical differencing on interpolated data):");
//        System.out.println("The refractive index of fused quartz at " + x1*1.0e9 + " nm is "+ y2);
//
//        double y3=cs.interpolate(x1);
//        System.out.println("Cubic spline:");
//        System.out.println("The refractive index of fused quartz at " + x1*1.0e9 + " nm is "+ y3);
//        System.out.println();
//
//        // Second interpolation at a wavelength of 590 nm
//        x1 = 5.9e-7;
//
//        y1=ci1.interpolate(x1);
//        System.out.println("Cubic interpolation (numerical differencing on entered data):");
//        System.out.println("The refractive index of fused quartz at " + x1*1.0e9 + " nm is "+ y1);
//
//        y2=ci2.interpolate(x1);
//        System.out.println("Cubic interpolation (numerical differencing on interpolated data):");
//        System.out.println("The refractive index of fused quartz at " + x1*1.0e9 + " nm is "+ y2);
//
//        y3=cs.interpolate(x1);
//        System.out.println("Cubic spline:");
//        System.out.println("The refractive index of fused quartz at " + x1*1.0e9 + " nm is "+ y3);
//        
//        x1 = 3.8e-7;
//
//        y1=ci1.interpolate(x1);
//        System.out.println("Cubic interpolation (numerical differencing on entered data):");
//        System.out.println("The refractive index of fused quartz at " + x1*1.0e9 + " nm is "+ y1);
//
//        y2=ci2.interpolate(x1);
//        System.out.println("Cubic interpolation (numerical differencing on interpolated data):");
//        System.out.println("The refractive index of fused quartz at " + x1*1.0e9 + " nm is "+ y2);
//
//        y3=cs.interpolate(x1);
//        System.out.println("Cubic spline:");
//        System.out.println("The refractive index of fused quartz at " + x1*1.0e9 + " nm is "+ y3);
//        
//        x1 = 4.6e-7;
//
//        y1=ci1.interpolate(x1);
//        System.out.println("Cubic interpolation (numerical differencing on entered data):");
//        System.out.println("The refractive index of fused quartz at " + x1*1.0e9 + " nm is "+ y1);
//
//        y2=ci2.interpolate(x1);
//        System.out.println("Cubic interpolation (numerical differencing on interpolated data):");
//        System.out.println("The refractive index of fused quartz at " + x1*1.0e9 + " nm is "+ y2);
//
//        y3=cs.interpolate(x1);
//        System.out.println("Cubic spline:");
//        System.out.println("The refractive index of fused quartz at " + x1*1.0e9 + " nm is "+ y3);
		
//		double[] x = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
//		double[] y = { 0.0, 0.8414, 0.9092, 0.1411, -0.7568, -0.9589, -0.2794, 0.6569, 0.9893, 0.4121, -0.5440 };
//		double[] xx = { 0, 0.2500, 0.500, 0.7500, 1, 1.2500, 1.500, 1.7500, 2, 2.2500, };
//		CubicInterpolation ci2 =new CubicInterpolation(x, y, 0);
//		for (int i = 0; i < xx.length; i++) {
//			System.out.println(ci2.interpolate(xx[i]));
//		}


	}
	
	
	public static double poly_interpolate(double[] dataX, double[] dataY, double x, int power)
	{
		int xIndex = 0;
		while (xIndex < dataX.length - (1 + power + (dataX.length - 1) % power) && dataX[xIndex + power] < x)
		{
			xIndex += power;
		}
 
		double matrix[][] = new double[power + 1][power + 2];
		for (int i = 0; i < power + 1; ++i)
		{
			for (int j = 0; j < power; ++j)
			{
				matrix[i][j] = Math.pow(dataX[xIndex + i], (power - j));
			}
			matrix[i][power] = 1;
			matrix[i][power + 1] = dataY[xIndex + i];
		}
		double[] coefficients = lin_solve(matrix);
		double answer = 0;
		for (int i = 0; i < coefficients.length; ++i)
		{
			answer += coefficients[i] * Math.pow(x, (power - i));
		}
		return answer;
	}
	
	
	// the gaussian elimnation algorithm
	// I won't explain this here because it's not in the scope of this tip
		private static double[] lin_solve(double[][] matrix)
		{
			double[] results = new double[matrix.length];
			int[] order = new int[matrix.length];
			for (int i = 0; i < order.length; ++i)
			{
				order[i] = i;
			}
			for (int i = 0; i < matrix.length; ++i)
			{
				// partial pivot
				int maxIndex = i;
				for (int j = i + 1; j < matrix.length; ++j)
				{
					if (Math.abs(matrix[maxIndex][i]) < Math.abs(matrix[j][i]))
					{
						maxIndex = j;
					}
				}
				if (maxIndex != i)
				{
					// swap order
					{
						int temp = order[i];
						order[i] = order[maxIndex];
						order[maxIndex] = temp;
					}
					// swap matrix
					for (int j = 0; j < matrix[0].length; ++j)
					{
						double temp = matrix[i][j];
						matrix[i][j] = matrix[maxIndex][j];
						matrix[maxIndex][j] = temp;
					}
				}
				if (Math.abs(matrix[i][i]) < 1e-15)
				{
					throw new RuntimeException("Singularity detected");
				}
				for (int j = i + 1; j < matrix.length; ++j)
				{
					double factor = matrix[j][i] / matrix[i][i];
					for (int k = i; k < matrix[0].length; ++k)
					{
						matrix[j][k] -= matrix[i][k] * factor;
					}
				}
			}
			for (int i = matrix.length - 1; i >= 0; --i)
			{
				// back substitute
				results[i] = matrix[i][matrix.length];
				for (int j = i + 1; j < matrix.length; ++j)
				{
					results[i] -= results[j] * matrix[i][j];
				}
				results[i] /= matrix[i][i];
			}
			double[] correctResults = new double[results.length];
			for (int i = 0; i < order.length; ++i)
			{
				// switch the order around back to the original order
				correctResults[order[i]] = results[i];
			}
			return results;
		}



}
