package org.jlab.trains;

// // CLAS Physics Imports
// import org.jlab.jnp.physics.*;

/**
 * 
 * Simple particle class with index added.
 *
 * @author mfmce
 * @version 15 Jul. 2021
 */

public class Part {

    private int idx, pid;
    private double px, py, pz;

    public Part(int idx, double px, double py, double pz) {
        this.idx = idx; this.px = px; this.py = py; this.pz = pz;
    }

    public Part(int pid, int idx, double px, double py, double pz) {
        this.pid = pid; this.idx = idx; this.px = px; this.py = py; this.pz = pz;
    }

    public int idx() { return this.idx; }
    public int pid() { return this.pid; }
    public double px() { return this.px; }
    public double py() { return this.py; }
    public double pz() { return this.pz; }
    public double pt() { return Math.sqrt(this.px*this.px+this.py*this.py); }
    public double p2() { return this.px*this.px+this.py*this.py+this.pz*this.pz; }
    public double p() { return Math.sqrt(this.px*this.px+this.py*this.py+this.pz*this.pz); }

} // class
