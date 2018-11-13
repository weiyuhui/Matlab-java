package CubicSpline;
/**
 * Representation of a cubic polynomial. 
 */ 
public final class Cubic { 
     
    // ******************************************************************** // 
    // Constructor. 
    // ******************************************************************** // 
 
    /**
     * Create a cubic polynomial of form a + b*x + c*x^2 + d*x^3. 
     *  
     * @param   a       A coefficient. 
     * @param   b       B coefficient. 
     * @param   c       C coefficient. 
     * @param   d       D coefficient. 
     */ 
    public Cubic(double a, double b, double c, double d) { 
        this.a = a; 
        this.b = b; 
        this.c = c; 
        this.d = d; 
    } 
 
     
    // ******************************************************************** // 
    // Evaluation. 
    // ******************************************************************** // 
 
    /**
     * Evaluate the polynomial for a given value. 
     *  
     * @param   x       X value to evaluate for. 
     * @return          The value of the polynomial for the given X. 
     */ 
    public double eval(double x) { 
        return ((d * x + c) * x + b) * x + a; 
    } 
 
     
    // ******************************************************************** // 
    // Private Data. 
    // ******************************************************************** // 
    
    // The coefficients. 
    private final double a, b, c, d; 
 
}

