package filtfilt;

import java.util.ArrayList;

public class ArrayList2Array {
	
	public static double[] toArray(ArrayList<Double> b) {
		double[] ret = new double[b.size()];
		for (int i = 0; i < ret.length; i++) {
			ret[i] = b.get(i);
		}
		return ret;
	}
	
	public static ArrayList<Double> toArrayList(double[] b) {
		ArrayList<Double> ret = new ArrayList<Double>();
		for (int i = 0; i < b.length; i++) {
			ret.add(b[i]);
		}
		return ret;
	}

}
