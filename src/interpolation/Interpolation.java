package interpolation;

import java.util.Arrays;
 
public class Interpolation {

    private int n;
    private Double[] xs;
    private Double[] ys;
     
    private boolean sp_initialized;
    private double[] sp_y2s;
 
    public Interpolation(Double[] _xs, Double[] _ys) {
        
        this.n = _xs.length;
        this.xs = Arrays.copyOf(_xs, _xs.length);
        this.ys = Arrays.copyOf(_ys, _ys.length);
       
        this.sp_initialized = false;
    } 
   
    public double spline(double x)
    {
        if (!this.sp_initialized) {
            // Assume Natural Spline Interpolation
            double p, qn, sig, un;
            double[] us;
            
            us = new double[n-1];
            sp_y2s = new double[n];
            us[0] = sp_y2s[0] = 0.0;
            
            for (int i=1; i<=n-2; i++) {
                sig = (xs[i] - xs[i-1]) / (xs[i+1] - xs[i-1]);
                p = sig * sp_y2s[i-1] + 2.0;
                sp_y2s[i] = (sig - 1.0) / p;
               us[i] = (ys[i+1] - ys[i]) / (xs[i+1] - xs[i]) - (ys[i] - ys[i-1]) / (xs[i] - xs[i-1]);
                us[i] = (6.0 * us[i] / (xs[i+1] - xs[i-1]) - sig * us[i-1]) / p;
            }
            qn = un = 0.0;
            
            sp_y2s[n-1] = (un - qn * us[n-2]) / (qn * sp_y2s[n-2] + 1.0);
            for (int k=n-2; k>=0; k--) {
                sp_y2s[k] = sp_y2s[k] * sp_y2s[k+1] + us[k];
            }

            this.sp_initialized = true;
        }
        
        int klo, khi, k;
        double h, b, a;
        
        klo = 0;
        khi = n-1;
        while (khi-klo > 1) {
            k = (khi+klo) >> 1;
            if (xs[k] > x)
                khi = k;
            else
                klo = k;
        }
        h = xs[khi] - xs[klo];
        if (h == 0.0) {
            throw new ArithmeticException();
        }
        a = (xs[khi] - x) / h;
        b = (x - xs[klo]) / h;
        return a*ys[klo] + b*ys[khi] + ((a*a*a-a)*sp_y2s[klo]+(b*b*b-b)*sp_y2s[khi])*(h*h)/6.0;
    }

 }