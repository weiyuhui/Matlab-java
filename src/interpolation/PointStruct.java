package interpolation;

public class PointStruct {
     
     public double dx;
     public double dy;
     
     public int ix;
     public int iy;
     
     public PointStruct(double dx, double dy) {
         this.dx = dx;
         this.dy = dy;
     }
 
     public PointStruct(int ix, int iy) {
         this.ix = ix;
         this.iy = iy;
     }
     
     public PointStruct(double dx, double dy,boolean round) {
         this.ix = RoundF(dx);
         this.iy = RoundF(dy);
     }
     
     public int RoundF(double a){
         return (int) Math.round(a);
     }
 }
