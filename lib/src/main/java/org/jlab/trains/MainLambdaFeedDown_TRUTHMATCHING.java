
// Java Imports
import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;

// CLAS Physics Imports
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriter;
import org.jlab.jnp.physics.*;

// J2ROOT Imports
import org.jlab.jroot.ROOTFile;
import org.jlab.jroot.TNtuple;

// CLAS QADB Import
import clasqa.QADB;

/**
* Test timing on bank reading in java.
*
* @version 19 July 2022
* @author Matthew McEneaney
*/

public class MainLambdaFeedDown_TRUTHMATCHING {

    //------------------------------ Constants ------------------------------//
    protected final static double beam_energy   = 10.6;
    protected final static double beam_mass     = 0.000510998950;
    protected final static double target_energy = 0.0;
    protected final static double target_mass   = 0.93827208816;

    // Constants
	private final static double mass_e      = (double) 0.0005110;		// [GeV/c^2]
	private final static double mass_p      = (double) 0.9383;			// [GeV/c^2]
	private final static double mass_pi     = (double) 0.13957;			// [GeV/c^2]
    private final static double pz_beam     = (double) 10.6;			// [GeV/c] RGA Fall 2018 Outbending beam energy
    private final static double mmin        = (double) 0.0;			    // [GeV/c^2] From Crystal Ball Fit
    private final static double mmax        = (double) 2.00;			// [GeV/c^2] From Crystal Ball Fit
    private final static double Q2min       = (double) 1.00;			// [-]
    private final static double Wmin        = (double) 2.00;			// [-]
    private final static double ymax        = (double) 0.80;		   	// [-]
    private final static double xFmin       = (double) -1.00;			// [-]
    private final static double zmax        = (double) 1.00;		   	// [-]

    // // built in lambdas
    // protected Cut _m  = (double m)  -> { return (m>mmin && m<mmax); };
    // protected Cut _Q2 = (double Q2) -> { return (Q2>Q2min); };
    // protected Cut _W  = (double W)  -> { return (W>Wmin);   };
    // protected Cut _y  = (double y)  -> { return (y<ymax);   };
    // protected Cut _xF = (double xF) -> { return (xF>xFmin); };
    // protected Cut _z  = (double z)  -> { return (z<zmax);   };

    protected int _counter   = 0;
    protected int _notify   = 10000;


	//------------------------------- Methods -------------------------------//
	public MainLambdaFeedDown_TRUTHMATCHING() {}

    protected Particle getParticle(Bank bank, int row) {
        int    pid = bank.getInt("pid",row);
        double px  = bank.getFloat("px",row);
        double py  = bank.getFloat("py",row);
        double pz  = bank.getFloat("pz",row);
        double vx  = bank.getFloat("vx",row);
        double vy  = bank.getFloat("vy",row);
        double vz  = bank.getFloat("vz",row);
        return new Particle(pid,px,py,pz,vx,vy,vz);
    } // protected Particle get Particle()

    protected Particle getElectron(Bank bank) {

        Particle electron = new Particle();

        // Loop bank
        for (int row=0; row<bank.getRows(); row++) {

            // Determine if particle is an electron and trigger particle and in FT or CD and chi2 cut
            int pid        = bank.getInt("pid",row);
            int status     = bank.getInt("status",row);
            double chi2pid = bank.getFloat("chi2pid",row);
        
            // if (pid == 11 && status<0 && Math.abs(chi2pid)<3 && Math.abs(status)>=2000 && Math.abs(status)<4000) {
            if (pid==11 && Math.abs(chi2pid)<3 && status<=-2000 /* status<=2000 *//*NOTE: NOT -2000???*/ /*&& status>-4000*/) {

                // Get electron
                Particle e_new = this.getParticle(bank,row);

                 // Reset maximum
                if (e_new.p() > electron.p()) {
                    electron = e_new;
                    if (row!=0) System.out.println("DEBUGGING: rowe!=0 instead "+row);//DEBUGGING
                }
            } // pid / status / chi2pid cuts
		} // for (int row=0; row<bank.getRows(); row++) {

        return electron;

    } // protected Particle getElectron()

    protected Particle getMCElectron(Bank bank) {

        // Get electron
        int row = 3;
        Particle electron = this.getParticle(bank,row);

        return electron;

    } // protected Particle getMCElectron()

    protected HashMap<String,Double> calculateKinematics(double helicity, Particle e, Particle p1, Particle p2, double chi2pide, double chi2pid1, double chi2pid2, double statuse, double status1, double status2) {

        // Set initial state lorentz vectors
        LorentzVector lv_beam   = new LorentzVector(); lv_beam.setPxPyPzM(0, 0, Math.sqrt(Math.pow(this.beam_energy,2) - Math.pow(this.beam_mass,2)), this.beam_mass); // Assumes all energy is along pz...?
        LorentzVector lv_target = new LorentzVector(); lv_target.setPxPyPzM(0, 0, 0, this.target_mass);

        // Set final state lorentz vectors
        LorentzVector lv_e      = e.vector();
        LorentzVector lv_q      = new LorentzVector(lv_beam); lv_q.sub(lv_e);
        LorentzVector lv_p1     = p1.vector();
        LorentzVector lv_p2     = p2.vector();
        LorentzVector lv_p      = new LorentzVector(lv_p1); lv_p.add(lv_p2);
        
        // Set gamma* nucleon frame lorentz vectors
        LorentzVector lv_gN        = new LorentzVector(lv_q);      lv_gN.add(lv_target);
        Vector3       gNBoost      = lv_gN.boostVector();          gNBoost.negative();
        LorentzVector lv_gN_target = new LorentzVector(lv_target); lv_gN_target.boost(gNBoost);
        LorentzVector lv_gN_e      = new LorentzVector(lv_e);      lv_gN_e.boost(gNBoost);
        LorentzVector lv_gN_q      = new LorentzVector(lv_q);      lv_gN_q.boost(gNBoost);
        LorentzVector lv_gN_p1     = new LorentzVector(lv_p1);     lv_gN_p1.boost(gNBoost);
        LorentzVector lv_gN_p2     = new LorentzVector(lv_p2);     lv_gN_p2.boost(gNBoost);
        LorentzVector lv_gN_p      = new LorentzVector(lv_p);      lv_gN_p.boost(gNBoost);

        // Set missing mass lorentz vectors
        LorentzVector lv_Mx    = new LorentzVector(lv_gN); lv_Mx.sub(lv_p);
        LorentzVector lv_Mx_p1 = new LorentzVector(lv_gN); lv_Mx_p1.sub(lv_p1);
        LorentzVector lv_Mx_p2 = new LorentzVector(lv_gN); lv_Mx_p2.sub(lv_p2);

        // Compute SIDIS variables
        double Q2   = -lv_q.mass2();
        double nu   = lv_q.e();
        double y    = lv_q.e() / lv_beam.e();
        double x    = Q2 / (2 * this.target_mass * nu);
        double W    = Math.sqrt(this.target_mass*this.target_mass+Q2 * (1 - x) / x);
        double Mh   = lv_p.mass();
        double Mx   = lv_Mx.mass();

        // Compute individual particle kinematics 1
        double z1      = lv_p1.e() / nu;
        double xF1     = lv_gN_p1.pz() / (W/2);
        double y1      = 1/2*Math.log((lv_p1.e()+lv_p1.pz())/(lv_p1.e()-lv_p1.pz()));
        double zeta1   = lv_gN_p1.e() / lv_gN_target.e();
        double Mx1     = lv_Mx_p1.mass();
        double phperp1 = lv_p1.vect().mag()*(lv_q.vect().cross(lv_p1.vect()).mag()/(lv_q.vect().mag()*lv_p1.vect().mag()));

        // Compute individual particle kinematics 2
        double z2      = lv_p2.e() / nu;
        double xF2     = lv_gN_p2.pz() / (W/2);
        double y2      = 1/2*Math.log((lv_p2.e()+lv_p2.pz())/(lv_p2.e()-lv_p2.pz()));
        double zeta2   = lv_gN_p2.e() / lv_gN_target.e();
        double Mx2     = lv_Mx_p2.mass();
        double phperp2 = lv_p2.vect().mag()*(lv_q.vect().cross(lv_p2.vect()).mag()/(lv_q.vect().mag()*lv_p2.vect().mag()));

        // Compute parent particle kinematics
        double zp      = lv_p.e() / nu;
        double xFp     = lv_gN_p.pz() / (W/2);
        double yp      = 1/2*Math.log((lv_p.e()+lv_p.pz())/(lv_p.e()-lv_p.pz()));
        double zetap   = lv_gN_p.e() / lv_gN_target.e();
        // double Mxp     = lv_Mx.mass();
        double phperpp = lv_p.vect().mag()*(lv_q.vect().cross(lv_p.vect()).mag()/(lv_q.vect().mag()*lv_p.vect().mag()));

        // Get phi_h_ in gN frame of momentum perpendicular to q relative to electron scattering plane ( just angle between q x l and q x p_hadron planes)
        Vector3 nhat    = lv_gN_q.vect().cross(lv_gN_e.vect());  // vA x vB
        Vector3 phihat  = lv_gN_q.vect().cross(lv_gN_p.vect()); // vC x vD
        double sign     = nhat.dot(lv_p.vect())>=0 ? 1 : -1;    // sign of (vA x vB) . vD
        double phi_h_p  = sign * Math.acos(nhat.dot(phihat)/(nhat.mag()*phihat.mag()));
        if (phi_h_p<0) phi_h_p = 2*Math.PI + phi_h_p;
        Vector3 phihat1 = lv_gN_q.vect().cross(lv_gN_p1.vect()); // vC x vD
        double sign1    = nhat.dot(lv_p1.vect())>=0 ? 1 : -1;    // sign of (vA x vB) . vD
        double phi_h_p1 = sign1 * Math.acos(nhat.dot(phihat1)/(nhat.mag()*phihat1.mag()));
        if (phi_h_p1<0) phi_h_p1 = 2*Math.PI + phi_h_p1;
        Vector3 phihat2 = lv_gN_q.vect().cross(lv_gN_p2.vect()); // vC x vD
        double sign2    = nhat.dot(lv_p2.vect())>=0 ? 1 : -1;    // sign of (vA x vB) . vD
        double phi_h_p2 = sign2 * Math.acos(nhat.dot(phihat2)/(nhat.mag()*phihat2.mag()));
        if (phi_h_p2<0) phi_h_p2 = 2*Math.PI + phi_h_p2;

        HashMap<String,Double> map = new HashMap<String,Double>();
        map.put("helicity",helicity);
        map.put("Q2",Q2);
        map.put("nu",nu);
        map.put("y",y);
        map.put("x",x);
        map.put("W",W);
        map.put("Mx",Mx);
        map.put("Mh",Mh);
        map.put("phi_h",phi_h_p);

        map.put("xF_p",xF1);
        map.put("y_p",y1);
        map.put("z_p",z1);
        map.put("zeta_p",zeta1);
        map.put("phperp_p",phperp1);
        map.put("phi_h_p",phi_h_p1);
        map.put("Mx_p",Mx1);

        map.put("xF_pim",xF2);
        map.put("y_pim",y2);
        map.put("z_pim",z2);
        map.put("zeta_pim",zeta2);
        map.put("phperp_pim",phperp2);
        map.put("phi_h_pim",phi_h_p2);
        map.put("Mx_pim",Mx2);

        map.put("xF_pa",xFp);
        map.put("y_pa",yp);
        map.put("z_pa",zp);
        map.put("zeta_pa",zetap);
        map.put("phperp_pa",phperpp);

        map.put("px_e",e.px());
        map.put("py_e",e.py());
        map.put("pz_e",e.pz());
        map.put("vx_e",e.vx());
        map.put("vy_e",e.vy());
        map.put("vz_e",e.vz());
        map.put("theta_e",e.theta());
        map.put("phi_e",e.phi());
        map.put("chi2pid_e",chi2pide);
        map.put("status_e",statuse);

        map.put("px_p",p1.px());
        map.put("py_p",p1.py());
        map.put("pz_p",p1.pz());
        map.put("vx_p",p1.vx());
        map.put("vy_p",p1.vy());
        map.put("vz_p",p1.vz());
        map.put("theta_p",p1.theta());
        map.put("phi_p",p1.phi());
        map.put("pid_p",(double)p1.pid());
        map.put("chi2pid_p",chi2pid1);
        map.put("status_p",status1);

        map.put("px_pim",p2.px());
        map.put("py_pim",p2.py());
        map.put("pz_pim",p2.pz());
        map.put("vx_pim",p2.vx());
        map.put("vy_pim",p2.vy());
        map.put("vz_pim",p2.vz());
        map.put("theta_pim",p2.theta());
        map.put("phi_pim",p2.phi());
        map.put("pid_pim",(double)p2.pid());
        map.put("chi2pid_pim",chi2pid2);
        map.put("status_pim",status2);

        map.put("px_pa",(double)lv_p.px());
        map.put("py_pa",(double)lv_p.py());
        map.put("pz_pa",(double)lv_p.pz());
        map.put("theta_pa",(double)lv_p.theta());
        map.put("phi_pa",(double)lv_p.phi());

        //----- BEGIN OLD -----//

        // double[] data = {
        //     helicity,
        //     Q2, nu, y, x, W, Mx, Mh, phi_h_p,
        //     xF1, y1, z1, zeta1, phperp1, phi_h_p1, Mx1, 
        //     xF2, y2, z2, zeta2, phperp2, phi_h_p2, Mx2,
        //     xFp, yp, zp, zetap, phperpp,
        //     e.px(), e.py(), e.pz(), e.vx(), e.vy(), e.vz(), e.theta(), e.phi(), chi2pide, statuse, 
        //     p1.px(), p1.py(), p1.pz(), p1.vx(), p1.vy(), p1.vz(), p1.theta(), p1.phi(), (double)p1.pid(), chi2pid1, status1, 
        //     p2.px(), p2.py(), p2.pz(), p2.vx(), p2.vy(), p2.vz(), p2.theta(), p2.phi(), (double)p2.pid(), chi2pid2, status2,
        //     lv_p.px(), lv_p.py(), lv_p.pz(), lv_p.theta(), lv_p.phi()
        // };

        // return data;
        //----- END OLD -----//

        return map;
    } // protected double[] calculateKinematics()

    protected HashMap<String,Double> calculateMCKinematics(double helicity, Particle e, Particle p1, Particle p2, Particle pa, double _ppid, double idx_pa, Particle ppa, double _pppid, double idx_ppa) {

        // Set initial state lorentz vectors
        LorentzVector lv_beam   = new LorentzVector(); lv_beam.setPxPyPzM(0, 0, Math.sqrt(Math.pow(this.beam_energy,2) - Math.pow(this.beam_mass,2)), this.beam_mass); // Assumes all energy is along pz...?
        LorentzVector lv_target = new LorentzVector(); lv_target.setPxPyPzM(0, 0, 0, this.target_mass);

        // Set final state lorentz vectors
        LorentzVector lv_e      = e.vector();
        LorentzVector lv_q      = new LorentzVector(lv_beam); lv_q.sub(lv_e);
        LorentzVector lv_p1     = p1.vector();
        LorentzVector lv_p2     = p2.vector();
        LorentzVector lv_p      = pa.vector();//new LorentzVector(lv_p1); lv_p.add(lv_p2);
        
        // Set gamma* nucleon frame lorentz vectors
        LorentzVector lv_gN        = new LorentzVector(lv_q);      lv_gN.add(lv_target);
        Vector3       gNBoost      = lv_gN.boostVector();          gNBoost.negative();
        LorentzVector lv_gN_target = new LorentzVector(lv_target); lv_gN_target.boost(gNBoost);
        LorentzVector lv_gN_e      = new LorentzVector(lv_e);      lv_gN_e.boost(gNBoost);
        LorentzVector lv_gN_q      = new LorentzVector(lv_q);      lv_gN_q.boost(gNBoost);
        LorentzVector lv_gN_p1     = new LorentzVector(lv_p1);     lv_gN_p1.boost(gNBoost);
        LorentzVector lv_gN_p2     = new LorentzVector(lv_p2);     lv_gN_p2.boost(gNBoost);
        LorentzVector lv_gN_p      = new LorentzVector(lv_p);      lv_gN_p.boost(gNBoost);

        // Set missing mass lorentz vectors
        LorentzVector lv_Mx    = new LorentzVector(lv_gN); lv_Mx.sub(lv_p);
        LorentzVector lv_Mx_p1 = new LorentzVector(lv_gN); lv_Mx_p1.sub(lv_p1);
        LorentzVector lv_Mx_p2 = new LorentzVector(lv_gN); lv_Mx_p2.sub(lv_p2);

        // Compute SIDIS variables
        double Q2   = -lv_q.mass2();
        double nu   = lv_q.e();
        double y    = lv_q.e() / lv_beam.e();
        double x    = Q2 / (2 * this.target_mass * nu);
        double W    = Math.sqrt(this.target_mass*this.target_mass+Q2 * (1 - x) / x);
        double Mh   = lv_p.mass();
        double Mx   = lv_Mx.mass();

        // Compute individual particle kinematics 1
        double z1      = lv_p1.e() / nu;
        double xF1     = lv_gN_p1.pz() / (W/2);
        double y1      = 1/2*Math.log((lv_p1.e()+lv_p1.pz())/(lv_p1.e()-lv_p1.pz()));
        double zeta1   = lv_gN_p1.e() / lv_gN_target.e();
        double Mx1     = lv_Mx_p1.mass();
        double phperp1 = lv_p1.vect().mag()*(lv_q.vect().cross(lv_p1.vect()).mag()/(lv_q.vect().mag()*lv_p1.vect().mag()));

        // Compute individual particle kinematics 2
        double z2      = lv_p2.e() / nu;
        double xF2     = lv_gN_p2.pz() / (W/2);
        double y2      = 1/2*Math.log((lv_p2.e()+lv_p2.pz())/(lv_p2.e()-lv_p2.pz()));
        double zeta2   = lv_gN_p2.e() / lv_gN_target.e();
        double Mx2     = lv_Mx_p2.mass();
        double phperp2 = lv_p2.vect().mag()*(lv_q.vect().cross(lv_p2.vect()).mag()/(lv_q.vect().mag()*lv_p2.vect().mag()));

        // Compute parent particle kinematics
        double zp      = lv_p.e() / nu;
        double xFp     = lv_gN_p.pz() / (W/2);
        double yp      = 1/2*Math.log((lv_p.e()+lv_p.pz())/(lv_p.e()-lv_p.pz()));
        double zetap   = lv_gN_p.e() / lv_gN_target.e();
        // double Mxp     = lv_Mx.mass();
        double phperpp = lv_p.vect().mag()*(lv_q.vect().cross(lv_p.vect()).mag()/(lv_q.vect().mag()*lv_p.vect().mag()));

        // Get phi_h_ in gN frame of momentum perpendicular to q relative to electron scattering plane ( just angle between q x l and q x p_hadron planes)
        Vector3 nhat    = lv_gN_q.vect().cross(lv_gN_e.vect());  // vA x vB
        Vector3 phihat  = lv_gN_q.vect().cross(lv_gN_p.vect()); // vC x vD
        double sign     = nhat.dot(lv_p.vect())>=0 ? 1 : -1;    // sign of (vA x vB) . vD
        double phi_h_p  = sign * Math.acos(nhat.dot(phihat)/(nhat.mag()*phihat.mag()));
        if (phi_h_p<0) phi_h_p = 2*Math.PI + phi_h_p;
        Vector3 phihat1 = lv_gN_q.vect().cross(lv_gN_p1.vect()); // vC x vD
        double sign1    = nhat.dot(lv_p1.vect())>=0 ? 1 : -1;    // sign of (vA x vB) . vD
        double phi_h_p1 = sign1 * Math.acos(nhat.dot(phihat1)/(nhat.mag()*phihat1.mag()));
        if (phi_h_p1<0) phi_h_p1 = 2*Math.PI + phi_h_p1;
        Vector3 phihat2 = lv_gN_q.vect().cross(lv_gN_p2.vect()); // vC x vD
        double sign2    = nhat.dot(lv_p2.vect())>=0 ? 1 : -1;    // sign of (vA x vB) . vD
        double phi_h_p2 = sign2 * Math.acos(nhat.dot(phihat2)/(nhat.mag()*phihat2.mag()));
        if (phi_h_p2<0) phi_h_p2 = 2*Math.PI + phi_h_p2;

        //TODO: json file, drop some entries...
        //TODO: get writeEvent method working...
        //TODO: rename MC entries all with _MC and idx_ppppid etc to be more obvious DONE
        //TODO: fix writeEvent method? should be close to fine...

        HashMap<String,Double> map = new HashMap<String,Double>();
        map.put("Q2_MC",Q2);
        map.put("nu_MC",nu);
        map.put("y_MC",y);
        map.put("x_MC",x);
        map.put("W_MC",W);
        map.put("Mx_MC",Mx);
        map.put("Mh_MC",Mh);
        map.put("phi_h_MC",phi_h_p);

        map.put("xF_p_MC",xF1);
        map.put("y_p_MC",y1);
        map.put("z_p_MC",z1);
        map.put("zeta_p_MC",zeta1);
        map.put("phperp_p_MC",phperp1);
        map.put("phi_h_p_MC",phi_h_p);
        map.put("Mx_p_MC",Mx1);

        map.put("xF_pim_MC",xF2);
        map.put("y_pim_MC",y2);
        map.put("z_pim_MC",z2);
        map.put("zeta_pim_MC",zeta2);
        map.put("phperp_pim_MC",phperp2);
        map.put("phi_h_pim_MC",phi_h_p2);
        map.put("Mx_pim_MC",Mx2);

        map.put("px_e_MC",e.px());
        map.put("py_e_MC",e.py());
        map.put("pz_e_MC",e.pz());
        map.put("vx_e_MC",e.vx());
        map.put("vy_e_MC",e.vy());
        map.put("vz_e_MC",e.vz());
        map.put("theta_e_MC",e.theta());
        map.put("phi_e_MC",e.phi());

        map.put("px_p_MC",p1.px());
        map.put("py_p_MC",p1.py());
        map.put("pz_p_MC",p1.pz());
        map.put("vx_p_MC",p1.vx());
        map.put("vy_p_MC",p1.vy());
        map.put("vz_p_MC",p1.vz());
        map.put("theta_p_MC",p1.theta());
        map.put("phi_p_MC",p1.phi());
        map.put("pid_p_MC",(double)p1.pid());

        map.put("px_pim_MC",p2.px());
        map.put("py_pim_MC",p2.py());
        map.put("pz_pim_MC",p2.pz());
        map.put("vx_pim_MC",p2.vx());
        map.put("vy_pim_MC",p2.vy());
        map.put("vz_pim_MC",p2.vz());
        map.put("theta_pim_MC",p2.theta());
        map.put("phi_pim_MC",p2.phi());
        map.put("pid_pim_MC",(double)p2.pid());

        map.put("pid_pa_MC",(double)_ppid);
        map.put("idx_pa_MC",(double)idx_pa);
        map.put("pid_ppa_MC",(double)_pppid);
        map.put("idx_ppa_MC",(double)idx_ppa);

        //----- BEGIN OLD -----//
        // Q2, nu, y, x, W, Mx, Mh, phi_h_p,
        // xF1, y1, z1, zeta1, phperp1, phi_h_p1, Mx1,
        // xF2, y2, z2, zeta2, phperp2, phi_h_p2, Mx2,
        // xFp, yp, zp, zetap, phperpp,
        // e.px(), e.py(), e.pz(), e.vx(), e.vy(), e.vz(), e.theta(), e.phi(),
        // p1.px(), p1.py(), p1.pz(), p1.vx(), p1.vy(), p1.vz(), p1.theta(), p1.phi(), (double)p1.pid(),
        // p2.px(), p2.py(), p2.pz(), p2.vx(), p2.vy(), p2.vz(), p2.theta(), p2.phi(), (double)p2.pid(),
        // pa.px(), pa.py(), pa.pz(), pa.vx(), pa.vy(), pa.vz(), pa.theta(), pa.phi(), (double)_ppid, idx_pa,
        // (double)_pppid, idx_ppa
        
        // return data;
        //----- END OLD -----//

        return map;
    } // protected double[] calculateMCKinematics()

	public static void main(String[] args) throws InterruptedException {

        if (args.length<1) { System.out.println("Usage: MainLambdaFeedDown_TRUTHMATCHING <inpath> [optional: <outpath>]");  return;}

        // PID search options
        boolean search_pid  = true;
        int     search_pid1 = 2212;
        int     search_pid2 = -211;

        // Charge search options
        boolean search_q    = false;
        int     search_q1   = 1;
        int     search_q2   = 1;

        // Read HIPO file
		String path = args[0];
		HipoReader reader = new HipoReader();
        reader.open(path);

        // Get HIPO REC::Particle bank
		Schema schema = reader.getSchemaFactory().getSchema("REC::Particle");
        Bank bank     = new Bank(schema);

        // Get HIPO MC::Lund bank
		Schema schema_mc = reader.getSchemaFactory().getSchema("MC::Lund");
        Bank bank_mc     = new Bank(schema_mc);

        // Get HIPO REC::Event bank
        Schema schema_event = reader.getSchemaFactory().getSchema("REC::Event");
        Bank bank_event     = new Bank(schema_event);

        // Create HIPO event for reading banks
        Event event   = new Event();

        // Open ROOT file for writing
        String outpath    = (String)(args.length>1 ? args[1] : "test.hipo");
        // ROOTFile rootfile = new ROOTFile(outpath);
        // String   tree     = "t";
        // String   title    = "title";
        // TNtuple  tntuple  = rootfile.makeNtuple(tree,title,
        //     "helicity:" +
        //     "Q2:nu:y:x:W:Mx:Mh:phi_h:" +
        //     "xF_p:y_p:z_p:zeta_p:phperp_p:phi_h_p:mx_p:" +
        //     "xF_pim:y_pim:z_pim:zeta_pim:phperp_pim:phi_h_pim:mx_pim:" +
        //     "xF_pa:y_pa:z_pa:zeta_pa:phperp_pa:" + //NOTE: ADDED
        //     "px_e:py_e:pz_e:vx_e:vy_e:vz_e:theta_e:phi_e:chi2pid_e:status_e:" +
        //     "px_p:py_p:pz_p:vx_p:vy_p:vz_p:theta_p:phi_p:pid_p:chi2pid_p:status_p:" +
        //     "px_pim:py_pim:pz_pim:vx_pim:vy_pim:vz_pim:theta_pim:phi_pim:pid_pim:chi2pid_pim:status_pim:"+
        //     "px_pa:py_pa:pz_pa:theta_pa:phi_pa:" + //NOTE: ADDED
        //     "Q2_MC:nu_MC:y_MC:x_MC:W_MC:Mx_MC:Mh_MC:phi_h_MC:" +
        //     "xF_p_MC:y_p_MC:z_p_MC:zeta_p_MC:phperp_p_MC:phi_h_p_MC:mx_p_MC:" +
        //     "xF_pim_MC:y_pim_MC:z_pim_MC:zeta_pim_MC:phperp_pim_MC:phi_h_pim_MC:mx_pim_MC:" +
        //     "xF_pa_MC:y_pa_MC:z_pa_MC:zeta_pa_MC:phperp_pa_MC:" + //NOTE: ADDED
        //     "px_e_MC:py_e_MC:pz_e_MC:vx_e_MC:vy_e_MC:vz_e_MC:theta_e_MC:phi_e_MC:" +
        //     "px_p_MC:py_p_MC:pz_p_MC:vx_p_MC:vy_p_MC:vz_p_MC:theta_p_MC:phi_p_MC:pid_p_MC:" +
        //     "px_pim_MC:py_pim_MC:pz_pim_MC:vx_pim_MC:vy_pim_MC:vz_pim_MC:theta_pim_MC:phi_pim_MC:pid_pim_MC:" +
        //     "px_pa_MC:py_pa_MC:pz_pa_MC:vx_pa_MC:vy_pa_MC:vz_pa_MC:theta_pa_MC:phi_pa_MC:pid_pa_MC:idx_pa_MC:pid_ppa_MC:idx_ppa_MC"); //NOTE: ADDED

        MainLambdaFeedDown_TRUTHMATCHING MainLambdaFeedDown_TRUTHMATCHING = new MainLambdaFeedDown_TRUTHMATCHING();

        // Open output HIPO file
        SchemaFactory schemaFactory = reader.getSchemaFactory(); // NOTE: For writing all event banks.

        schemaFactory.readFile("/Users/mfm45/ppim_kinbank.json"); //NOTE: Need absolute path here.
        HipoWriter writer = new HipoWriter(schemaFactory);
        System.out.println("keys = "+schemaFactory.getSchemaKeys()); //NOTE: For just writing specific banks.
        schemaFactory.getSchema("MC::TruthMatching").show();
        writer.open(outpath);

        
        //----- BEGIN DEBUGGING -----//
        // boolean __rrg__ = false; //while (__rrg__) { //
        // int __max_events__ = 50000000;
        int __counter__ = 0;
        int __write_counter__ = 0;
        //----- END DEBUGGING -----//

        // Loop events
		while(reader.hasNext()==true) { //NOTE: COMMENTED OUT FOR DEBUGGING

            reader.nextEvent(event);

            //----- BEGIN DEBUGGING -----//
            __counter__++;
            // if (__counter__>__max_events__) break;
            //----- END DEBUGGING -----//

            // Define data list of hashmaps to which to add bank data
            ArrayList<HashMap<String,Double>> data = new ArrayList<HashMap<String,Double>>();

            // Get event helicity
            event.read(bank_event);
            double helicity = (double)bank_event.getInt("helicity",0);

            // Get scattered electron REC
            event.read(bank); //NOTE: Read bank IMPORTANT!
            Particle electron = MainLambdaFeedDown_TRUTHMATCHING.getElectron(bank); //NOTE: Since we require status<0 rowe == 0 always
            if (electron.p()==0) continue;
            int rowe = 0;

            // Get scattered electron MC
            event.read(bank_mc); //NOTE: Read bank IMPORTANT!
            Particle _electron = MainLambdaFeedDown_TRUTHMATCHING.getMCElectron(bank_mc); //NOTE: Since we require status<0 rowe == 0 always
            if (_electron.p()==0) continue;
            int _rowe = 0;

            event.read(bank); //NOTE: Read REC::Particle bank IMPORTANT!
            for (int row1=0; row1<bank.getRows(); row1++) { // Row loop 1
                for (int row2=0; row2<bank.getRows(); row2++) { // Row loop 2
                    if (row2==row1 || row1==rowe || row2==rowe) continue; //NOTE: Uniqueness requirement
					double pid1  = bank.getInt("pid", row1);
                    double pid2  = bank.getInt("pid", row2);
                    if (search_pid && (pid1!=search_pid1 || pid2!=search_pid2)) continue; //NOTE: Check pids
                    double q1    = bank.getByte("charge", row1);
                    double q2    = bank.getByte("charge", row2);
                    if (search_q && (q1!=search_q1 || q2!=search_q2)) continue; //NOTE: Check charges
                    
                    // Get bank info for particles
                    Particle p1 = MainLambdaFeedDown_TRUTHMATCHING.getParticle(bank,row1);
                    Particle p2 = MainLambdaFeedDown_TRUTHMATCHING.getParticle(bank,row2);

                    // Also get chi2pid and status
                    double chi2pide = bank.getFloat("chi2pid",rowe);
                    double statuse  = (double)bank.getInt("status",rowe);
                    double chi2pid1 = bank.getFloat("chi2pid",row1);
                    double status1  = (double)bank.getInt("status",row1);
                    double chi2pid2 = bank.getFloat("chi2pid",row2);
                    double status2  = (double)bank.getInt("status",row2);

                    // Get kinematics
                    HashMap<String,Double> data_rec = MainLambdaFeedDown_TRUTHMATCHING.calculateKinematics(helicity,electron,p1,p2,chi2pide,chi2pid1,chi2pid2,statuse,status1,status2);

                    // Now Loop MC bank for matching decay
                    event.read(bank_mc); //NOTE: Read MC::Lund bank IMPORTANT!
                    for (int _row1=3; _row1<bank_mc.getRows(); _row1++) { // Row loop 1 //NOTE: BEGIN AT 3 FOR MC SINCE BEAM TARGET VIRTUAL PHOTON AND SCATTERED ELECTRON ARE FIRST.
                        for (int _row2=3; _row2<bank_mc.getRows(); _row2++) { // Row loop 2
                            if (_row2==_row1 || _row1==_rowe || _row2==_rowe) continue; //NOTE: Uniqueness requirement
                            double _pid1  = bank_mc.getInt("pid", _row1);
                            double _pid2  = bank_mc.getInt("pid", _row2);
                            if (search_pid && (_pid1!=search_pid1 || _pid2!=search_pid2)) continue; //NOTE: Check pids
                            // double _q1    = bank_mc.getByte("charge", _row1);
                            // double _q2    = bank_mc.getByte("charge", _row2);
                            // if (search_q && (q1!=search_q1 || q2!=search_q2)) continue; //NOTE: Check charges
                            
                            // Get bank info for particles
                            Particle _p1 = MainLambdaFeedDown_TRUTHMATCHING.getParticle(bank_mc,_row1);
                            Particle _p2 = MainLambdaFeedDown_TRUTHMATCHING.getParticle(bank_mc,_row2);

                            // Also get parent indices
                            int _pidx1 = bank_mc.getInt("parent",_row1);
                            int _pidx2 = bank_mc.getInt("parent",_row2);
                            // if (_pidx1!=_pidx2) { continue; }//DEBUGGING: TODO //NOTE: COMMENTED OUT SINCE WE WANT ALL BACKGROUND!
                            double _ppid = bank_mc.getInt("pid",_pidx1-1); //NOTE: MC::Lund index begins at 1 but bank index begins at 0.
                            Particle parent = MainLambdaFeedDown_TRUTHMATCHING.getParticle(bank_mc,_pidx1-1);
                            int _ppidx = bank_mc.getInt("parent",_pidx1-1);
                            double _pppid = bank_mc.getInt("pid",_ppidx-1);
                            Particle pparent = MainLambdaFeedDown_TRUTHMATCHING.getParticle(bank_mc,_ppidx-1);

                            // Get MC Kinematics
                            HashMap<String,Double> data_mc = MainLambdaFeedDown_TRUTHMATCHING.calculateMCKinematics(helicity,_electron,_p1,_p2,parent,_ppid,_pidx1-1.0,pparent,_pppid,_ppidx-1.0);

                            // Add data to output file
                            HashMap<String,Double> data_new = new HashMap<String,Double>();
                            data_new.putAll(data_rec);
                            data_new.putAll(data_mc);//NOTE: THIS RELIES ON KEYS BEING UNIQUE!!!
                            data_new.put("idx_p",(double)row1);
                            data_new.put("idx_pim",(double)row2);
                            data_new.put("idx_p_MC",(double)_row1);
                            data_new.put("idx_pim_MC",(double)_row2);
                            data_new.put("idx_pa_pim_MC",(double)_pidx2);
                            data.add(data_new);//NOTE: This appends to a list of hashmaps that should be defined in the outermost loops scope.
                            
                            // MainLambdaFeedDown_TRUTHMATCHING.writeEvent(schemaFactory,writer,event,__row__,data_new);//TODO!!!

                        } // for (int row2=0; row2<bank.getRows(); row2++) { //NOTE: MC LOOP
                    } // for (int row1

                    event.read(bank);//NOTE: read bank IMPORTANT!
                
                } // for (int row2=0; row2<bank.getRows(); row2++) { //NOTE: REC LOOP
            } // for (int row1=0; row1<bank.getRows(); row1++) {

            // Write event if you added data
            MainLambdaFeedDown_TRUTHMATCHING.writeEvent(schemaFactory,bank,event,writer,data);
            __write_counter__++;

		} // while(reader.hasNext()==true)

        // // Close ROOT file
        // tntuple.write();
        // rootfile.close();

        // Close HIPO writer
        writer.close();

        System.out.println("DEBUGGING: __counter__ = "+__counter__);//DEBUGGING
        System.out.println("DEBUGGING: __write_counter__ = "+__write_counter__);//DEBUGGING
        
    } // Main()

    private void writeEvent(SchemaFactory factory, Bank bank, Event event, HipoWriter writer, ArrayList<HashMap<String,Double>> data) {
    
        Bank newbank = new Bank(factory.getSchema("MC::TruthMatching"),data.size()); //IMPORTANT: DO NOT LEAVE OFF SECOND ARGUMENT: data.size() = # of rows 
        int row = 0;
        boolean added_data = false;
        for (HashMap<String,Double> __data__ : data) {
            if (__data__.size()==0) continue;
                        if (__data__.get("Mh")>mmax || 
                __data__.get("Q2")<Q2min || 
                __data__.get("W")<Wmin || 
                __data__.get("y")>ymax || 
                __data__.get("xF_pa")<xFmin || 
                __data__.get("z_pa")>zmax ) continue;
            this.fillKinematics(newbank, row, __data__);
            added_data = true;
            // event.reset(); //NOTE: For just writing specific banks.
            row++;
        } // for (HashMap<String,Double> ...)
        if (data.size()>0 && added_data) {
            // event.write(bank); //TODO: Write list of banks...
            event.write(newbank);
            writer.addEvent(event);
        }
    } // writeEvent()

    // public HashMap<String,Boolean> computeCuts(HashMap<String,Double> kinematics) {

    //     HashMap<String,Boolean> cuts = new HashMap<String,Boolean>();

    //     // Process entries
    //     cuts.put("mass",this._m.cut(kinematics.get("mass")));
    //     cuts.put("Q2",this._Q2.cut(kinematics.get("Q2")));
    //     cuts.put("W",this._W.cut(kinematics.get("W")));
    //     cuts.put("y",this._y.cut(kinematics.get("y")));
    //     cuts.put("xF",this._xF.cut(kinematics.get("xF")));
    //     cuts.put("z",this._z.cut(kinematics.get("z")));

    //     return cuts;
    // } // computeCuts()

    public void fillKinematics(Bank bank, int row, HashMap<String,Double> kinematics) {

        for(String key : kinematics.keySet()) {
            if ((key.contains("pid") && !key.contains("chi2")) || key.contains("helicity") || key.contains("status") || key.contains("idx")) { bank.putInt(key,row,(int)kinematics.get(key).intValue()); }
            else                    { bank.putDouble(key,row,(double)kinematics.get(key)); }
        }
    } // fillKinematics()

} // class   


// System.out.println("keyset = ");
//             int counter__=0;
//             for (String key : __data__.keySet()) {
//                 System.out.print(" "+key+(counter__%5==0 ? " ,\n" : " , "));
//                 counter__++;
//             } System.out.println("}");
//             System.out.println("DEBUGGING: __data__.containsKey(\"Mh\") = "+__data__.containsKey("Mh"));//DEBUGGING
//             System.out.println("DEBUGGING: __data__.containsKey(\"Q2\") = "+__data__.containsKey("Q2"));//DEBUGGING
//             System.out.println("DEBUGGING: __data__.containsKey(\"W\") = "+__data__.containsKey("W"));//DEBUGGING
//             System.out.println("DEBUGGING: __data__.containsKey(\"y\") = "+__data__.containsKey("y"));//DEBUGGING
//             System.out.println("DEBUGGING: __data__.containsKey(\"xF_pa\") = "+__data__.containsKey("xF_pa"));//DEBUGGING
//             System.out.println("DEBUGGING: __data__.containsKey(\"z_pa\") = "+__data__.containsKey("z_pa"));//DEBUGGING

