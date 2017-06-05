package Orbital;

import java.util.ArrayList;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display.Circle;

public class Particle extends Circle {
	double distance;
	double velocityx;
	double velocityy;
	double accelerationx;
	double accelerationy;
	double mass;


	public double getXpos() {
		return x;
	}






	public void setXpos(double xpos) {
		this.x = xpos;
	}







	public double getYpos() {
		return y;
	}







	public void setYpos(double ypos) {
		this.y = ypos;
	}







	public double getDistance() {
		return distance;
	}







	public void setDistance(double distance) {
		this.distance = distance;
	}







	public double getVelocityx() {
		return velocityx;
	}







	public void setVelocityx(double velocityx) {
		this.velocityx = velocityx;
	}







	public double getVelocityy() {
		return velocityy;
	}







	public void setVelocityy(double velocityy) {
		this.velocityy = velocityy;
	}







	public double getAccelerationx() {
		return accelerationx;
	}







	public void setAccelerationx(double accelerationx) {
		this.accelerationx = accelerationx;
	}







	public double getAccelerationy() {
		return accelerationy;
	}







	public void setAccelerationy(double accelerationy) {
		this.accelerationy = accelerationy;
	}







	public double getMass() {
		return mass;
	}







	public void setMass(double mass) {
		this.mass = mass;
	}







	public void step(double gravityconstant, double timestep, ArrayList<Particle> particles) {
		double accelerationx = 0;
		double accelerationy = 0;
		for (int i = 0; i < particles.size(); i++) {
			if (particles.get(i) != this) {
				double massother = particles.get(i).getMass();
				double xposother = particles.get(i).getXpos();
				double yposother = particles.get(i).getYpos();
				double angle = Math.atan2(yposother-y,xposother-x);
				double distance = Math.sqrt((yposother-y)*(yposother-y) + (xposother-x)*(xposother-x));
				double force = (gravityconstant*mass*massother)/(Math.pow(distance, 2));

				accelerationx += (Math.cos(angle))* (force/mass);
				accelerationy += (Math.sin(angle))*(force/mass);
			}
		}
		setAccelerationx(accelerationx);
		setAccelerationy(accelerationy);
		setVelocityx(timestep* getAccelerationx() + getVelocityx());
		setVelocityy(timestep* getAccelerationy() + getVelocityy());
		setXpos(getXpos() + getVelocityx() *timestep);
		setYpos(getYpos() + getVelocityy()*timestep);
		setXY(getXpos(), getYpos());

	}
}
