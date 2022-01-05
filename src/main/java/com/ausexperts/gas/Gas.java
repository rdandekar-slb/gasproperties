package com.ausexperts.gas;

public class Gas {
    private double specificGravity=0.65;
    private double moleFractionCO2=0.0;
    private double moleFractionH2S=0.0;
    //private double moleFractionN2=0.0;


    public double gaszfactor=0.0;
    public double gasfvf=0.0;
    public double gasviscosity=0.0;
    public double gasdensity=0.0;

    

    public Gas(){    }

    public Gas(double spgr){
        this.specificGravity=spgr;
    }

    public Gas(double spgr, double co2, double h2s /*, double n2*/){
        this.specificGravity=spgr;
        this.moleFractionCO2=co2;
        this.moleFractionH2S=h2s;
        //this.moleFractionN2=n2;
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

    public double getSpecificGravity(){
        return this.specificGravity;
    }
    /*
    public void setMoleFractionN2(double moleFractionN2) {
        this.moleFractionN2 = moleFractionN2;
    }
    */
    
    private void getproperties(double pressure, double temperature) {
        double Tpc=168+325*this.specificGravity-12.5*this.specificGravity*this.specificGravity;
        double ppc=677+15.0*this.specificGravity-37.5*this.specificGravity*this.specificGravity;
        double epsilon=120*(Math.pow(this.moleFractionCO2+this.moleFractionH2S,0.9)-Math.pow(this.moleFractionH2S, 1.6))+15*(Math.pow(this.moleFractionH2S, 0.5)-Math.pow(this.moleFractionH2S, 4.0));
        Tpc=Tpc-epsilon;
        ppc=(ppc*Tpc)/(Tpc+this.moleFractionH2S*(1-this.moleFractionH2S)*epsilon);
        double Tpr=(temperature+460)/Tpc;
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
        double r3=a6 + a7/Tpr + a8/(Tpr*Tpr);
        double r4=a9*(a7/Tpr + a8/(Tpr*Tpr))  ;
        double r5=a10/(Tpr*Tpr*Tpr);
        double rhor=0;
        double rhor_new=0;
        double f_rhor=0;
        double fprime_rhor=0;
        
        for (int i=0;i<50;i++){
            if (i==0){
                rhor=0.27*ppr/Tpr;
            }
            f_rhor=r1*rhor - r2/rhor + r3*rhor*rhor - r4*Math.pow(rhor,5) + r5*(1+a11*rhor*rhor)*rhor*rhor*Math.exp(-a11*rhor*rhor) + 1;
            fprime_rhor=r1 + r2/(rhor*rhor) + 2*r3*rhor - 5*r4*Math.pow(rhor, 4) + 2*r5*rhor*Math.exp(-a11*rhor*rhor)*((1+2*a11*rhor*rhor*rhor) - a11*rhor*rhor*(1+a11*rhor*rhor));
            rhor_new=rhor - f_rhor/fprime_rhor;

            if (Math.abs(rhor - rhor_new) < 1e-12){
                rhor=rhor_new;
                break;
            }
            rhor=rhor_new;
        }
        double apparentMolecularWeight=28.96*this.specificGravity;
        this.gaszfactor=0.27*ppr/(rhor*Tpr);
        this.gasdensity=(apparentMolecularWeight*pressure) / (this.gaszfactor*10.73*(temperature+460));
        this.gasfvf=0.02827*this.gaszfactor*(temperature+460)/pressure;
        double k=(9.4+0.02*apparentMolecularWeight)*Math.pow(temperature+460, 1.5)/(209+19*apparentMolecularWeight+(temperature+460));
        double x=3.5+(986/(temperature+460))+0.01*apparentMolecularWeight;
        double y=2.4-0.2*x;
        this.gasviscosity=0.0001*k*Math.exp(x*Math.pow(this.gasdensity/62.4, y));
        return;
    }

    public static void main(String[] args) {
        Gas g = new Gas(0.72,0,0);
        g.getproperties(2000 ,  140);
        System.out.println(g.gaszfactor);
        System.out.println(g.gasfvf);
        System.out.println(g.gasviscosity);
        System.out.println(g.gasdensity);
        
    }



}
