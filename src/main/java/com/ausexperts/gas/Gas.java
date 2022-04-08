package com.ausexperts.gas;

public class Gas {
    private double specificGravity = 0.65;
    private double moleFractionCO2 = 0.0;
    private double moleFractionH2S = 0.0;
    // private double moleFractionN2=0.0;

    private double gaszfactor = 0.0;
    private double gasfvf = 0.0;
    private double gasviscosity = 0.0;

    public Gas() {
    }

    public Gas(double spgr) {
        this.specificGravity = spgr;
    }

    public Gas(double spgr, double co2, double h2s /* , double n2 */) {
        this.specificGravity = spgr;
        this.moleFractionCO2 = co2;
        this.moleFractionH2S = h2s;
        // this.moleFractionN2=n2;
    }

    public void setSpecificGravity(double specificGravity) {
        this.specificGravity = specificGravity;
    }

    public void setMoleFractionCO2(double moleFractionCO2) {
        this.moleFractionCO2 = moleFractionCO2;
    }

    public void setMoleFractionH2S(double moleFractionH2S) {
        this.moleFractionH2S = moleFractionH2S;
    }
    /*
     * public void setMoleFractionN2(double moleFractionN2) {
     * this.moleFractionN2 = moleFractionN2;
     * }
     */

    public double getGaszfactor(double pressure, double temperature) {
        double Tpc=168+325*specificGravity-12.5*specificGravity*specificGravity;
        double ppc=677+15.0*specificGravity-37.5*specificGravity*specificGravity;
        double epsilon=120*(Math.pow(moleFractionCO2+moleFractionH2S,0.9)-Math.pow(moleFractionH2S, 1.6))+15*(Math.pow(moleFractionH2S, 0.5)-Math.pow(moleFractionH2S, 4.0));
        Tpc-=epsilon;
        ppc=(ppc*Tpc)/(Tpc+moleFractionH2S*(1-moleFractionH2S)*epsilon);
        double Tpr=(temperature+459.67)/Tpc;
        double ppr=pressure/ppc;


        double a1=0.3265;
        double a2=-1.0700;
        double a3=-0.5339;
        double a4=0.01569;
        double a5=-0.05165;
        double a6=0.5475;
        double a7=-0.7361;
        double a8=0.1844;
        double a9=0.1056;
        double a10=0.6134;
        double a11=0.7210;

        double r1=a1 + a2/Tpr + a3/Math.pow(Tpr, 3) + a4/Math.pow(Tpr, 4) + a5/Math.pow(Tpr, 5);
        double r2=0.27*ppr/Tpr;
        double r3=a6 + a7/Tpr + a8/Math.pow(Tpr, 2);
        double r4=a9*(a7/Tpr + a8/Math.pow(Tpr, 2));
        double r5=a10/Math.pow(Tpr, 3);
        

        double rhor_k=0.27*ppr/Tpr;
        double gaszfactor=1;
        double rhor_k_1,f_rhor, f_prime_rhor;
        while (true){
        f_rhor= r1*rhor_k - r2/rhor_k + r3*Math.pow(rhor_k, 2) - r4*Math.pow(rhor_k, 5) + r5*(1+a11*Math.pow(rhor_k, 2))*Math.pow(rhor_k, 2)*Math.exp(-a11*rhor_k*rhor_k)+1;
        f_prime_rhor = r1 + r2/Math.pow(rhor_k,2) + 2*r3*rhor_k - 5*r4*Math.pow(rhor_k, 4) + 2*r5*rhor_k*Math.exp(-a11 *Math.pow(rhor_k, 2))*((1+2*a11*Math.pow(rhor_k, 3))- a11*Math.pow(rhor_k, 2)*(1+a11*Math.pow(rhor_k, 2)));
        rhor_k_1 = rhor_k - f_rhor/f_prime_rhor;
            if (Math.abs(rhor_k - rhor_k_1) < 1.0e-8){
                break;
                }
            else{
                rhor_k = rhor_k_1;
            }
        }
        gaszfactor = 0.27*ppr/(rhor_k*Tpr);

        return gaszfactor;
    }

}
