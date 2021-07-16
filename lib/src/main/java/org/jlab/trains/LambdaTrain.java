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
import org.jlab.jnp.hipo.data.HipoNode;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriter;
import org.jlab.jnp.physics.*;

/**
 * 
 * Standalone skim for Λ -> p+ + π- analysis. (Shortened version with invariant mass and kinematic cuts)
 *
 * @author mfmce
 * @version 15 Jul. 2021
 */

public class LambdaTrain {

    public LambdaTrain() {}

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

    // built in lambdas
    protected Cut _m  = (double m)  -> { return (m>mmin && m<mmax); };
    protected Cut _Q2 = (double Q2) -> { return (Q2>Q2min); };
    protected Cut _W  = (double W)  -> { return (W>Wmin);   };
    protected Cut _y  = (double y)  -> { return (y<ymax);   };
    protected Cut _xF = (double xF) -> { return (xF>xFmin); };
    protected Cut _z  = (double z)  -> { return (z<zmax);   };

    protected int _counter   = 0;
    protected int _notify   = 10000;

    public static void main(String[] args){

        if (args.length == 0) { System.out.println("Usage: LambdaTrain [infile] [outfile]"); return; }

        LambdaTrain test_wagon = new LambdaTrain();

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
            SchemaFactory schemaFactory = reader.getSchemaFactory(); // NOTE: For writing all event banks.
            // SchemaFactory schemaFactory = new SchemaFactory(); //NOTE: For just writing specific banks.
            // schemaFactory.addSchema(reader.getSchemaFactory().getSchema("REC::Particle")); //NOTE: For just writing specific banks.
            // schemaFactory.addSchema(reader.getSchemaFactory().getSchema("RUN::config")); //NOTE: For just writing specific banks.
            // schemaFactory.addSchema(reader.getSchemaFactory().getSchema("REC::Event")); //NOTE: For just writing specific banks.
            // schemaFactory.addSchema(reader.getSchemaFactory().getSchema("MC::Lund")); //NOTE: For just writing specific banks.
            schemaFactory.readFile("/CLAS12-Trains/lib/etc/kinematics.json"); //NOTE: Need absolute path here.
            // System.out.println("keys = "+schemaFactory.getSchemaKeys()); //NOTE: For just writing specific banks.
            HipoWriter writer = new HipoWriter(schemaFactory);
            writer.open(outpath);

            // Loop events
            while(reader.hasNext()==true) {
                reader.nextEvent(event);
                test_wagon.processDataEvent(event, schemaFactory, writer);
            }
            System.out.println("counter = "+test_wagon._counter);
            writer.close();
        }
    } // main()
	
	public void processDataEvent(Event event, SchemaFactory factory, HipoWriter writer) {

		Bank bank = new Bank(factory.getSchema("REC::Particle"));
		event.read(bank);

        // Initialize outer scope variables
		boolean flag_e  = false; boolean flag_p = false; boolean flag_pi = false;
        ArrayList<ArrayList<Part>> data = new ArrayList<ArrayList<Part>>();
        Part electron = new Part(0,0,0,0);

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
                    double px_e = bank.getFloat("px", index);
                    double py_e = bank.getFloat("py", index);
                    double pz_e = bank.getFloat("pz", index);
                    electron = new Part(index,px_e,py_e,pz_e);
                }
            } // if (pid == 11)
        } // for (index in bank...)

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
                        flag_p = true; flag_pi = true;

                        // Get momenta for proton
                        double px_p = bank.getFloat("px", index);
                        double py_p = bank.getFloat("py", index);
                        double pz_p = bank.getFloat("pz", index);
                        Part proton = new Part(index,px_p,py_p,pz_p);

                        // Get momenta for pion
                        double px_pi = bank.getFloat("px", index2);
                        double py_pi = bank.getFloat("py", index2);
                        double pz_pi = bank.getFloat("pz", index2);
                        Part pion = new Part(index2,px_pi,py_pi,pz_pi);

                        // Add data to arrays
                        ArrayList<Part> particles = new ArrayList<Part>();
                        particles.add(proton); particles.add(pion);
                        data.add(particles);
                        
                    } // if (pid2 == -211)
                } // for ( index2 in bank...)
            } // if (pid == 2212)
		} // for ( index in bank...)

        // Create new event to write to
        Event nevent = event;       //NOTE: For writing all banks.
        // Event nevent = new Event(); //NOTE: For just writing specific banks.

        if (flag_e && flag_p && flag_pi) { this.writeEvent(factory,bank,nevent,writer,electron,data); }
	} // processDataEvent()

    private void writeEvent(SchemaFactory factory, Bank bank, Event event, HipoWriter writer, Part electron, ArrayList<ArrayList<Part>> data) {
    
        ArrayList<HashMap<String,Double>> lkinematics = new ArrayList<HashMap<String,Double>>();
        for (ArrayList<Part> particles : data) {
            HashMap<String,Double> kinematics = this.computeKinematics(electron, particles.get(0), particles.get(1));
            HashMap<String,Boolean> cuts      = this.computeCuts(kinematics);

            if (cuts.get("mass") && cuts.get("Q2") && cuts.get("W") && cuts.get("y") && cuts.get("xF") && cuts.get("z")) {
                lkinematics.add(kinematics);
            }
        } // for (ArrayList<Double> ...)
        int row = 0;
        Bank kinBank = new Bank(factory.getSchema("REC::Kinematics"),lkinematics.size()); //IMPORTANT: DO NOT LEAVE OFF SECOND ARGUMENT: data.size() = # of rows 
        for (HashMap<String,Double> kinematics : lkinematics) {

            this.fillKinematics(kinBank, row, kinematics);
            // event.reset(); //NOTE: For just writing specific banks.
            if (row == 0) { this._counter++; }
            if ((this._counter % this._notify) == 0 && row == 0) { System.out.println("PROCESSED "+this._counter+" EVENTS"); }
            row++;
        } // for (HashMap<String,Double> ...)
        if (lkinematics.size()>0) {
            // event.write(bank); //TODO: Write list of banks...
            event.write(kinBank);
            writer.addEvent(event);
        }
    } // writeEvent()

    private HashMap<String,Double> computeKinematics(Part electron, Part proton, Part pion) {

        double px_e = electron.px(); double py_e = electron.py(); double pz_e = electron.pz(); double p2_e = electron.p2();
        double px_p = proton.px(); double py_p = proton.py(); double pz_p = proton.pz(); double p2_p = proton.p2();
        double px_pi = pion.px(); double py_pi = pion.py(); double pz_pi = pion.pz(); double p2_pi = pion.p2();
        
        // Initialize map
        HashMap<String, Double> kinematics = new HashMap<String,Double>();
        kinematics.put("idxe",(double)electron.idx());
        kinematics.put("idxp",(double)proton.idx());
        kinematics.put("idxpi",(double)pion.idx());

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

        kinematics.put("mass",m);
        kinematics.put("Q2",Q2);
        kinematics.put("nu",nu);
        kinematics.put("W",W);
        kinematics.put("x",x);
        kinematics.put("y",y);
        kinematics.put("xF",xF);
        kinematics.put("z",z);

        return kinematics;
    } // computeKinematics()

    public HashMap<String,Boolean> computeCuts(HashMap<String,Double> kinematics) {

        HashMap<String,Boolean> cuts = new HashMap<String,Boolean>();

        // Process entries
        cuts.put("mass",this._m.cut(kinematics.get("mass")));
        cuts.put("Q2",this._Q2.cut(kinematics.get("Q2")));
        cuts.put("W",this._W.cut(kinematics.get("W")));
        cuts.put("y",this._y.cut(kinematics.get("y")));
        cuts.put("xF",this._xF.cut(kinematics.get("xF")));
        cuts.put("z",this._z.cut(kinematics.get("z")));

        return cuts;
    } // computeCuts()

    public void fillKinematics(Bank bank, int row, HashMap<String,Double> kinematics) {

        for(String key : kinematics.keySet()) {
            if (key.contains("id")) { bank.putInt(key,row,kinematics.get(key).intValue()); }
            else                    { bank.putDouble(key,row,kinematics.get(key)); }
        }
    } // fillKinematics()

} // class
