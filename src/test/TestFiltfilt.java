package test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import biz.source_code.dsp.filter.FilterPassType;
import biz.source_code.dsp.filter.IirFilterCoefficients;
import biz.source_code.dsp.filter.IirFilterDesignExstrom;
import filtfilt.ArrayList2Array;
import filtfilt.FilterFilterArray;
import filtfilt.Filtfilt;
import filtfilt.Mat;

public class TestFiltfilt {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		// Filtfilt.doFiltfilt();
//		int[] vec = {1,2,3};
//		for (Integer integer : vec) {
//			System.out.println(integer);
//		}
		
		int order = 5;
		double Wn1 = 0.15;
		double Wn2 = 0.5;
		IirFilterCoefficients iirFilterCoefficients;        
		iirFilterCoefficients = IirFilterDesignExstrom.design(FilterPassType.lowpass, order, Wn1/2,Wn2/2);
		double[] a = iirFilterCoefficients.a;
		double[] b = iirFilterCoefficients.b;
		
		ArrayList<Double> B = new ArrayList<Double>();
		ArrayList<Double> A = new ArrayList<Double>();
		
		B = ArrayList2Array.toArrayList(b);
		A = ArrayList2Array.toArrayList(a);
		
		
		ArrayList<Double> X = new ArrayList<Double>();
		
//		B.add(4.0);
//		B.add(5.0);
//		B.add(6.0);
//		
//		A.add(1.0);
//		A.add(2.0);
//		A.add(3.0);
		
		X.add(1.0);
		X.add(2.0);
		X.add(3.0);
		X.add(4.0);
		X.add(5.0);
		X.add(6.0);
		X.add(7.0);
		X.add(8.0);
		X.add(9.0);
		X.add(10.0);
		X.add(11.0);
		X.add(12.0);
		X.add(13.0);
		X.add(14.0);
		X.add(15.0);
		X.add(16.0);
		X.add(17.0);
		X.add(18.0);
		
		
		ArrayList<Double> y = Filtfilt.doFiltfilt(B,A,X);
		for (int i = 0; i < y.size(); i++)
			System.out.println(y.get(i));
		
		
//		double[] y2  = doFiltfilt(B1, A1, X1);
//		
//		for (int i = 0; i < y2.length; i++)
//			System.out.println(y2[i]);
		
	}
	
	


	

}
