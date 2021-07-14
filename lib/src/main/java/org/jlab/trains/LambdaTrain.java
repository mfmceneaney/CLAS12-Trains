package org.jlab.trains;

// Java Imports
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

// CLAS Physics Imports
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriter;
import org.jlab.jnp.physics.*;

/**
 * 
 * Standalone skim for Λ -> p+ + π- analysis. (Shortened version with invariant mass and kinematic cuts)
 *
 * @author mfmce
 * @version 14 Jul. 2021
 */

public class LambdaTrain {

    public LambdaTrain() {}

    // Constants
	private final static double mass_e      = (double) 0.0005110;		// [GeV/c^2]
	private final static double mass_p      = (double) 0.9383;			// [GeV/c^2]
	private final static double mass_pi     = (double) 0.13957;			// [GeV/c^2]
    private final static double pz_beam     = (double) 10.6;			// [GeV/c] RGA Fall 2018 Outbending beam energy
    private final static double mmin        = (double) 1.10;			// [GeV/c^2] From Crystal Ball Fit
    private final static double mmax        = (double) 1.13;			// [GeV/c^2] From Crystal Ball Fit
    private final static double Q2min       = (double) 1.00;			// [-]
    private final static double Wmin        = (double) 2.00;			// [-]
    private final static double ymax        = (double) 0.80;		   	// [-]
    private final static double xFmin       = (double) -1.00;			// [-]
    private final static double zmax        = (double) 1.00;		   	// [-]

    // built in lambdas
    protected Cut _m  = (double m)  -> { return (m>mmin && m<mmax); };
    protected Cut _Q2 = (double Q2) -> { return (Q2>Q2min); };
    protected Cut _W  = (double W)  -> { return (W>Wmin);   };
    protected Cut _y  = (double y)  -> { return (y<ymax);   };
    protected Cut _xF = (double xF) -> { return (xF>xFmin); };
    protected Cut _z  = (double z)  -> { return (z<zmax);   };

    protected int counter   = 0;
    protected int _notify   = 100;

    public static void main(String[] args){

        if (args.length == 0) { System.out.println("Usage: Train [infile] [outfile]"); return; }

        Train test_wagon = new Train();

        // Check input file/directory for HIPO files
        String path = new String(args[0]);
        String outpath = new String(args[1]);
        File folder = new File(path);

        // Single file
        if (folder.isFile() && folder.getName().contains(".hipo")) {

            // Open input HIPO file
            HipoReader reader = new HipoReader();
            reader.open(path);
            Event event   = new Event();

            // Open output HIPO file
            HipoWriter writer = new HipoWriter(reader.getSchemaFactory());
            writer.open(outpath);

            while(reader.hasNext()==true) {

                reader.nextEvent(event);
                boolean flag = test_wagon.processDataEvent(event, reader.getSchemaFactory(), writer);
            }
            System.out.println("counter = "+test_wagon.counter);

            writer.close();
        }
    } // main()
	
	public boolean processDataEvent(Event event, SchemaFactory factory, HipoWriter writer) {

		Bank bank = new Bank(factory.getSchema("REC::Particle"));
		event.read(bank);

		boolean flag_e  = false;
		boolean flag_pi = false;
        boolean flag_p  = false;
        boolean flag_m  = false;
        boolean flag_Q2 = false;
        boolean flag_W  = false;
        boolean flag_y  = false;
        boolean flag_xF = false;
        boolean flag_z  = false;
        double px_e, py_e, pz_e, p2_e, px_pi, py_pi, pz_pi, p2_pi, px_p, py_p, pz_p, p2_p;
        px_e = 0; py_e = 0; pz_e = 0; p2_e = 0;
        px_pi = 0; py_pi = 0; pz_pi = 0; p2_pi = 0;
        px_p = 0; py_p = 0; pz_p = 0; p2_p = 0;

        for (int index = 0; index < bank.getRows(); index++) {
			int pid = bank.getInt("pid", index);
			
			if (pid == 11) {
                int status = bank.getInt("status", index);
                double chi2pid = (double) bank.getFloat("chi2pid", index);
                boolean safety = flag_e;
				flag_e = (status>-3000 && status<=-1000 && Math.abs(chi2pid)<3);
                if (safety && flag_e) System.out.println(" *** WARNING *** Double scattered electron event!");
                if (flag_e) {
                    // Get momenta for scattered electron
                    px_e = bank.getFloat("px", index);
                    py_e = bank.getFloat("py", index);
                    pz_e = bank.getFloat("pz", index);
                    p2_e = (px_e*px_e+py_e*py_e+pz_e*pz_e);
                }
            }
        }

        // Loop for protons
		for (int index = 0; index < bank.getRows(); index++) {
            if (!flag_e) continue;
			int pid = bank.getInt("pid", index);
            if (pid == 2212) {

                // Loop for pions
                for (int index2 = 0; index2 < bank.getRows(); index2++) {
                    int pid2 = bank.getInt("pid", index2);
                    if (pid2 == -211) { 

                        // Set flags
                        flag_p  = true; flag_pi = true;

                        // Get momenta for proton
                        px_p = bank.getFloat("px", index);
                        py_p = bank.getFloat("py", index);
                        pz_p = bank.getFloat("pz", index);
                        p2_p = (px_p*px_p+py_p*py_p+pz_p*pz_p);

                        // Get momenta for pion
                        px_pi = bank.getFloat("px", index2);
                        py_pi = bank.getFloat("py", index2);
                        pz_pi = bank.getFloat("pz", index2);
                        p2_pi = (px_pi*px_pi+py_pi*py_pi+pz_pi*pz_pi);  

                        HashMap<String, Boolean> cuts = this.computeCuts(px_e,py_e,pz_e,p2_e,px_pi,py_pi,pz_pi,p2_pi,px_p,py_p,pz_p,p2_p);
                        flag_m  = cuts.get("m");
                        flag_Q2 = cuts.get("Q2");
                        flag_W  = cuts.get("W");
                        flag_y  = cuts.get("y");
                        flag_xF = cuts.get("xF");
                        flag_z  = cuts.get("z");

                        if (flag_m && flag_Q2 && flag_W && flag_y && flag_xF && flag_z) {
                            counter++;
                            if ((counter % this._notify) == 0) { System.out.println("PROCESSED "+counter+" EVENTS"); }
                            writer.addEvent(event);
                            return true;
                        }
                    } // if (pid2 == -211)
                } // for ( index2 in bank...)   
            } // if (pid == 2212)        
		} // for ( index in bank...)

		return false;
	} // processDataEvent()

    private HashMap<String, Boolean> computeCuts(double px_e, double py_e, double pz_e, double p2_e,
                                                      double px_pi, double py_pi, double pz_pi, double p2_pi,
                                                      double px_p, double py_p, double pz_p, double p2_p) {
        HashMap<String, Boolean> cuts = new HashMap<String, Boolean>();

        // Set interaction lorentz vectors
        LorentzVector lv_max = new LorentzVector();
        LorentzVector lv_beam = new LorentzVector();
        LorentzVector lv_target = new LorentzVector();
        LorentzVector lv_L0 = new LorentzVector();

        lv_max.setPxPyPzM(px_e, py_e, pz_e, mass_e);		// Set vectors in Lab Frame
        lv_beam.setPxPyPzM(0.0, 0.0, pz_beam, mass_e);	    // Incoming beam at 10.6 GeV for RGA Fall 2018
        lv_target.setPxPyPzM(0.0, 0.0, 0.0, mass_p);
        lv_L0.setPxPyPzE(px_p+px_pi,py_p+py_pi,pz_p+pz_pi,
            Math.sqrt(p2_p+mass_p*mass_p)+Math.sqrt(p2_pi+mass_pi*mass_pi)); // Fixed from m+ + m- to E+ + E- 12/17/20

        LorentzVector q = new LorentzVector(lv_beam);
        q.sub(lv_max);
        LorentzVector gN = new LorentzVector(q);
        gN.add(lv_target);
        Vector3 gNBoost = gN.boostVector();
        gNBoost.negative();
        LorentzVector boostedLambda = new LorentzVector(lv_L0);
        boostedLambda.boost(gNBoost);
        LorentzVector boostedMax = new LorentzVector(lv_max);
        boostedMax.boost(gNBoost);

        // Compute SIDIS variables
        double Q2 = (-1) * (q.mass2());
        double nu = q.e();
        double z  = lv_L0.e() / nu;
        double y  = q.e() / lv_beam.e();		
        double x  = Q2 / (2 * mass_p * nu);
        double W  = Math.sqrt(mass_p*mass_p+Q2 * (1 - x) / x);
        double xF = boostedLambda.pz() / (W/2);

        // Energy and mass
        double E_p  = Math.sqrt(p2_p+Math.pow(mass_p,2));
        double E_pi = Math.sqrt(p2_pi+Math.pow(mass_pi,2));
        double mom2 = Math.pow(px_p+px_pi,2) + Math.pow(py_p+py_pi,2) + Math.pow(pz_p+pz_pi,2);
        double m    = Math.sqrt(Math.pow(E_p+E_pi,2)-mom2);

        cuts.put("m",this._m.cut(m));
        cuts.put("Q2",this._Q2.cut(Q2));
        cuts.put("W",this._W.cut(W));
        cuts.put("y",this._y.cut(y));
        cuts.put("xF",this._xF.cut(xF));
        cuts.put("z",this._z.cut(z));

        return cuts;
    } // computeCuts()

} // class
