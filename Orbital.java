package Orbital;


import java.awt.Color;
import java.util.ArrayList;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display.*;
import org.opensourcephysics.frames.*;

//import Orbital.Particle;

public class Orbital extends AbstractSimulation {
	ArrayList<Trail> trails = new ArrayList<Trail>(); //creates array list for trails
	ArrayList<Particle> particles = new ArrayList<Particle>(); // creates array list for particles
	double maxX = Double.NEGATIVE_INFINITY; //sets double maxX to negative infinity for proof of Kepler's 1st Law, Kepler's 2nd Law, and Bonus 
	double minX = Double.POSITIVE_INFINITY; //sets double minX to positive infinity for proof of Kepler's 1st Law  Kepler's 2nd Law, and Bonus
	double maxY = Double.NEGATIVE_INFINITY;  //sets double maxY to negative infinity for proof of Kepler's 1st Law, Kepler's 2nd Law, and Bonus
	double minY = Double.POSITIVE_INFINITY; //sets double minX to negative infinity for proof of Kepler's 1st Law, Kepler's 2nd Law, and Bonus
	double maxsX = Double.NEGATIVE_INFINITY; //sets double maxsX to negative infinity for proof of Bonus
	double minsX = Double.POSITIVE_INFINITY; //sets double minsX to positive infinity for proof of Bonus
	double maxsY = Double.NEGATIVE_INFINITY; //sets double maxsY to negative infinity for proof of Bonus
	double minsY = Double.POSITIVE_INFINITY; //sets double minsY to negative infinity for proof of Bonus
	ArrayList<Double> xhistory = new ArrayList<Double>(); //creates array list xhistory for proof of Kepler's 2nd Law, Kepler's 3rd Law, and Bonus
	ArrayList<Double> yhistory = new ArrayList<Double>(); //creates array list yhistory for proof of Kepler's 2nd Law, Kepler's 3rd Law, and Bonus
	ArrayList<Double> xhistorysun = new ArrayList<Double>(); //creates array list xhistorysun for proof of Bonus
	ArrayList<Double> yhistorysun = new ArrayList<Double>(); //creates arraylist yhistorysun for proof of Bonus
	Trail triangle = new Trail(); //create new trail triangle for proof of Kepler's 2nd Law
	int cycles = 0; //declare int cycles, set it to zero
	int cyclessun = 0; //declare int cyclessun, set it to zero
	int pcycles = 0; //declare int pcycles, set it to zero
	double period = 0; //declare double period, set it to zero
	double gravityconstant = 6.67E-11; //declare double gravityconstant, set it to zero
	Particle sun = new Particle(); //create new particle sun
	Particle earth = new Particle(); //create new particle earth
	int step = 0; //declare int step, set it to zero
	int sunstep = 0;
	DisplayFrame frame = new DisplayFrame("X" , "Y", "Display Frame Test"); //create frame

	/**
	 * This method calls on the step method in the Particle Class for each particle. Proves Kepler's First Law, Kepler's Second Law, Kepler's Third Law, and Bonus. 
	 */
	protected void doStep() {
		step +=1; //increment step by one
		sunstep+=1;
		earth.step(gravityconstant, control.getDouble("timestep"), particles); //calls on step method in the Particle Class for particle earth with parameters gravityconstant, timestep, and particles
		sun.step(gravityconstant, control.getDouble("timestep"), particles); //calls on step method in the Particle Class for particle sun with parameters gravityconstant, timestep, and particles

		Trail suntrail = trails.get(0); //set suntrail to the first trail in trails array list
		Trail earthtrail = trails.get(1); //set earthtrail to the second trail in trails array list
		suntrail.addPoint(sun.getX(), sun.getY()); //add point to suntrail
		earthtrail.addPoint(earth.getX(), earth.getY()); //add point to earthtrail

		if (control.getBoolean("Kepler's First Law") == true) { //if "Kepler's First Law" is set to true in control panel
			if(particles.get(1).getXpos() > maxX){ //if the x position of the particle is greater than the recorded maximum x position
				maxX = particles.get(1).getXpos(); 
			}
			if(particles.get(1).getXpos() < minX){
				minX = particles.get(1).getXpos(); 
			}
			if(particles.get(1).getYpos() > maxY){
				maxY = particles.get(1).getYpos(); 
			}
			if(particles.get(1).getYpos() < minY){
				minY = particles.get(1).getYpos();
			}

			xhistory.add(particles.get(1).getXpos());
			yhistory.add(particles.get(1).getYpos());
			if(step > 10) {
				if(particles.get(1).getYpos() > yhistory.get(xhistory.size()-10) && particles.get(1).getXpos() > xhistory.get(xhistory.size()-10)) {
					cycles +=1;
				}
			}

			if (cycles >= 1) {

				double h = (maxX + minX)/2;
				double k = (maxY + minY)/2;
				double a = (maxX - minX)/2;
				double b = (maxY - minY)/2;
				double xpos = particles.get(1).getXpos();
				double ypos = particles.get(1).getYpos();

				Trail ellipse = new Trail();
				ellipse.color = Color.red;
				frame.addDrawable(ellipse);

				if(a>b) {
					if (ypos > 0) {
						ellipse.addPoint(xpos, k + Math.sqrt(b*b - (b*b*(xpos-h)*(xpos-h))/(a*a)));
					}
					if (ypos < 0) {
						ellipse.addPoint(xpos, k - Math.sqrt(b*b - (b*b*(xpos-h)*(xpos-h))/(a*a)));
					}
				}
				if(a<b) {
					if (ypos > 0) {
						ellipse.addPoint(xpos, k + Math.sqrt(a*a - (a*a*(xpos-h)*(xpos-h))/(b*b)));
					}
					if (ypos < 0) {
						ellipse.addPoint(xpos, k - Math.sqrt(a*a - (a*a*(xpos-h)*(xpos-h))/(b*b)));
					}
				}

			}

		}

		//PROVE KEPLER'S SECOND LAW: Is the areal velocity constant?
		if (control.getBoolean("Kepler's Second Law") == true) { //if "Kepler's Second Law" is set to true in control panel
			xhistory.add(particles.get(1).getXpos());
			yhistory.add(particles.get(1).getYpos());

			frame.removeDrawable(triangle);
			triangle.clear();
			PlotFrame timevsarealvelocity = new PlotFrame("time","arealvelocity", "time vs. areal velocity");
			timevsarealvelocity.setVisible(true);

			if(step > 10) {
				triangle.addPoint(xhistory.get(xhistory.size()-10), yhistory.get(yhistory.size()-10));
				triangle.addPoint(particles.get(1).getXpos(), particles.get(1).getYpos());
				triangle.addPoint(particles.get(0).getXpos(), particles.get(0).getYpos());

				frame.addDrawable(triangle);

				double x1 = xhistory.get(xhistory.size()-10);
				double y1 = yhistory.get(yhistory.size()-10);
				double x2 = particles.get(1).getXpos();
				double y2 = particles.get(1).getYpos();
				double x3 = particles.get(0).getXpos();
				double y3 = particles.get(0).getYpos();

				double trianglearea = Math.abs((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2);
				System.out.println(trianglearea);

			}

		}

		//PROVE KEPLER'S THIRD LAW: How is the period related to the major axis?
		if (control.getBoolean("Kepler's Third Law") == true) { //if "Kepler's Third Law" is set to true in control panel
			if(particles.get(1).getXpos() > maxX){
				maxX = particles.get(1).getXpos(); 
			}
			if(particles.get(1).getXpos() < minX){
				minX = particles.get(1).getXpos(); 
			}
			if(particles.get(1).getYpos() > maxY){
				maxY = particles.get(1).getYpos(); 
			}
			if(particles.get(1).getYpos() < minY){
				minY = particles.get(1).getYpos();
			}
			xhistory.add(particles.get(1).getXpos());
			yhistory.add(particles.get(1).getYpos());
			if(step > 10) {
				if(particles.get(1).getYpos() > yhistory.get(yhistory.size()-10) && particles.get(1).getXpos() > xhistory.get(xhistory.size()-10)) {
					cycles +=1;
				}
			}

			double a;
			double b;
			double major;
			xhistory.add(particles.get(1).getXpos());
			yhistory.add(particles.get(1).getYpos());

			if (step > 10) {
				if(yhistory.get(yhistory.size()-10) < 0 && particles.get(1).getYpos() > yhistory.get(yhistory.size()-10) && particles.get(1).getXpos() < xhistory.get(xhistory.size()-10)) {
					pcycles +=1;
				}
			}
			if(pcycles < 1) {
				period = period + control.getDouble("timestep");
			}
			else {
				System.out.print("period : semimajor axis = " + period + " : ");
				if(cycles>=1) {
					a = (maxX - minX)/2;
					b = (maxY - minY)/2;
					if(a>b) {
						major = a;
						System.out.println(major);
					}
					else {
						major = b;
						System.out.print(major);
					}
					System.out.println("(GMsunT^2)/(4π^2a^3) = " + (gravityconstant*period*period*sun.getMass())/(4*Math.PI*Math.PI*major*major*major));

				}
			}
		}

		//BONUS: If the mass of the planet is approaching the mass of the sun, will the planet make an ellipse around the star? Will the star make an ellipse?
		if (control.getBoolean("Bonus") == true) { //if "Bonus" is set to true in control panel
			if(particles.get(1).getXpos() > maxX){
				maxX = particles.get(1).getXpos(); 
			}
			if(particles.get(1).getXpos() < minX){
				minX = particles.get(1).getXpos(); 
			}
			if(particles.get(1).getYpos() > maxY){
				maxY = particles.get(1).getYpos(); 
			}
			if(particles.get(1).getYpos() < minY){
				minY = particles.get(1).getYpos();
			}

			xhistory.add(particles.get(1).getXpos());
			yhistory.add(particles.get(1).getYpos());

			if(step > 10) {
				//System.out.println(particles.get(1).getYpos() + ">" + yhistory.get(xhistory.size()-10));
				if(particles.get(1).getYpos() > yhistory.get(yhistory.size()-10) && particles.get(1).getXpos() > xhistory.get(xhistory.size()-10)) {
					//cycles = true;
					cycles +=1;
					//System.out.println(cycles);
				}
			}

			//center of ellipse:
			if (cycles >= 1) {
				double h = (maxX + minX)/2;
				double k = (maxY + minY)/2;
				double a = (maxX - minX)/2;
				double b = (maxY - minY)/2;
				double c1;
				double c2;
				double xpos = particles.get(1).getXpos();
				double ypos = particles.get(1).getYpos();

				Trail ellipse = new Trail();
				ellipse.color = Color.red;
				frame.addDrawable(ellipse);
				Circle circle1 = new Circle();
				Circle circle2 = new Circle();
				frame.addDrawable(circle1);
				frame.addDrawable(circle2);
				if(a>b) {
					c1 = Math.sqrt(a*a-b*b);
					c2 = -Math.sqrt(a*a-b*b);
					circle1.setXY(h + c1, k);
					circle2.setXY(h + c2, k);
					System.out.println("earth's eccentricity:" + Math.sqrt(a*a + b*b)/a);
					if (ypos > 0) {
						ellipse.addPoint(xpos, k + Math.sqrt(b*b - (b*b*(xpos-h)*(xpos-h))/(a*a)));
					}
					if (ypos < 0) {
						ellipse.addPoint(xpos, k - Math.sqrt(b*b - (b*b*(xpos-h)*(xpos-h))/(a*a)));
					}
				}
				if(a<b) {
					c1 = Math.sqrt(b*b-a*a);
					c2 = -Math.sqrt(b*b-a*a);
					circle1.setXY(h + c1, k);
					circle2.setXY(h + c2, k);
					System.out.println("earth's eccentricity:" + Math.sqrt(a*a + b*b)/a);
					if (ypos > 0) {
						ellipse.addPoint(xpos, k + Math.sqrt(a*a - (a*a*(xpos-h)*(xpos-h))/(b*b)));
					}
					if (ypos < 0) {
						ellipse.addPoint(xpos, k - Math.sqrt(a*a - (a*a*(xpos-h)*(xpos-h))/(b*b)));
					}
				}

			}
			
			if(particles.get(0).getXpos() > maxsX){
				maxsX = particles.get(0).getXpos(); 
			}
			if(particles.get(0).getXpos() < minsX){
				minsX = particles.get(0).getXpos(); 
			}
			if(particles.get(0).getYpos() > maxsY){
				maxsY = particles.get(0).getYpos(); 
			}
			if(particles.get(0).getYpos() < minsY){
				minsY = particles.get(0).getYpos();
			}




			xhistorysun.add(particles.get(0).getXpos());
			yhistorysun.add(particles.get(0).getYpos());


			if(sunstep > 10) {
				//System.out.println(sunstep);
				//System.out.println("youre in");
				//System.out.println("compare" + particles.get(0).getYpos() + " // " + yhistorysun.get(yhistorysun.size()-1) + "//" + yhistorysun.get(yhistorysun.size()-2));
				if(particles.get(0).getYpos() < yhistorysun.get(yhistorysun.size()-2) && particles.get(0).getXpos() < xhistorysun.get(xhistorysun.size()-2)) {
					cyclessun += 1;
				}
			}



			if (cyclessun >= 1) {
				double hsun = (maxsX + minsX)/2;
				double ksun = (maxsY + minsY)/2;
				double asun = (maxsX - minsX)/2;
				double bsun = (maxsY - minsY)/2;
				double c1;
				double c2;
				double xpossun = particles.get(0).getXpos();
				double ypossun = particles.get(0).getYpos();

				Trail ellipsesun = new Trail();
				ellipsesun.color = Color.green;
				frame.addDrawable(ellipsesun);
				Circle circle1 = new Circle();
				Circle circle2 = new Circle();
				circle1.color = Color.MAGENTA;
				circle2.color = Color.MAGENTA;
				frame.addDrawable(circle1);
				frame.addDrawable(circle2);
				if(asun>bsun) {
					c1 = Math.sqrt(asun*asun-bsun*bsun);
					c2 = -Math.sqrt(asun*asun-bsun*bsun);
					circle1.setXY(hsun + c1, ksun);
					circle2.setXY(hsun + c2, ksun);
					System.out.println("sun's eccentricity:" + Math.sqrt(asun*asun + bsun*bsun)/asun);
					if (ypossun > 0) {
						ellipsesun.addPoint(xpossun, ksun + Math.sqrt(bsun*bsun - (bsun*bsun*(xpossun-hsun)*(xpossun-hsun))/(asun*asun)));
						
					}
					if (ypossun < 0) {
						ellipsesun.addPoint(xpossun, ksun - Math.sqrt(bsun*bsun - (bsun*bsun*(xpossun-hsun)*(xpossun-hsun))/(asun*asun)));
						
					}
				}
				if(asun<bsun) {
					c1 = Math.sqrt(asun*asun-bsun*bsun);
					c2 = -Math.sqrt(asun*asun-bsun*bsun);
					circle1.setXY(hsun + c1, ksun);
					circle2.setXY(hsun + c2, ksun);
					System.out.println("sun's eccentricity:" + Math.sqrt(asun*asun + bsun*bsun)/asun);
					if (ypossun > 0) {
						ellipsesun.addPoint(xpossun, ksun + Math.sqrt(asun*asun - (asun*asun*(xpossun-hsun)*(xpossun-hsun))/(bsun*bsun)));
						
					}
					if (ypossun < 0) {
						ellipsesun.addPoint(xpossun, ksun - Math.sqrt(asun*asun - (asun*asun*(xpossun-hsun)*(xpossun-hsun))/(bsun*bsun)));
						
					}
				}

			}

		}


	}

	public void initialize()
	{
		sun.setXpos(control.getDouble("x position sun"));
		sun.setYpos(control.getDouble("y position sun"));
		sun.setVelocityx(control.getDouble("x velocity sun"));
		sun.setVelocityy(control.getDouble("y velocity sun"));
		sun.setMass(control.getDouble("mass sun"));
		earth.setXpos(control.getDouble("x position earth")); 
		earth.setYpos(control.getDouble("y position earth"));
		earth.setVelocityx(control.getDouble("x velocity earth"));
		earth.setVelocityy(control.getDouble("y velocity earth"));
		earth.setMass(control.getDouble("mass earth"));

		earth.color = Color.blue;
		sun.color = Color.yellow;

		particles.add(sun);
		particles.add(earth);



		Trail suntrail = new Trail();
		Trail earthtrail = new Trail();
		trails.add(earthtrail);
		frame.addDrawable(earthtrail);
		trails.add(suntrail);
		frame.addDrawable(suntrail);
		frame.setVisible(true);
		frame.addDrawable(earth);
		frame.addDrawable(sun);

		frame.setPreferredMinMax(-2E11, 2E11, -2E11, 2E11);
	}

	public void reset() {
		control.setValue("x position sun", 0);
		control.setValue("y position sun", 0);
		control.setValue("x position earth", 1.5E11);
		control.setValue("y position earth", 0);
		control.setValue("mass sun", 1.989E30);
		control.setValue("mass earth", 5.972E24);
		control.setValue("x velocity sun", 0);
		control.setValue("y velocity sun", 0);
		control.setValue("x velocity earth", 0);
		control.setValue("y velocity earth", 20000); //from 3
		control.setValue("timestep", 1000); //make one less for second law
		control.setValue("Kepler's First Law", false);
		control.setValue("Kepler's Second Law", false);
		control.setValue("Kepler's Third Law", false);
		control.setValue("Bonus", false);
		//for bonus, mass earth 2E29, velocity sun is -2000
	}

	public static void main(String[] args) { //main
		SimulationControl.createApp(new Orbital()); //Creates a SIP animation control and establishes communication between the control and the model
	}
}
