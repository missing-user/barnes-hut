#pragma once

#define GLM_FORCE_PURE
//#define GLM_PRECISION_MEDIUMP_FLOAT
#include "ofMain.h"
#include "ofxGui.h"

#include "Simulation/Distributions.h"
#include "Simulation/Order.h"
#include "OctTree/Particle.h"
#include "Simulation/Simulation.h"

class ofApp : public ofBaseApp{
	private:
		std::vector<Particle> particles;
  	void initializeParticles();
		void drawBoxes();
	public:
		void setup();
		void update();
		void draw();
		
		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y);
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
		
		ofMesh mesh;
		ofIcoSpherePrimitive icoSphere;
		ofEasyCam cam;
		
  ofxPanel gui;

		

  // Adjustable parameters
  ofParameter<double> timestep_slider;
  ofParameter<int> max_depth_slider;
  ofParameter<int> max_per_node_slider;
  ofParameter<int> num_particles_slider;
  ofParameter<myfloat> theta_slider;
  ofParameter<myfloat> mass_slider;
  ofParameter<bool> brute_force_toggle;

  // stats
  ofParameter<int> min_depth_slider;
  ofParameter<bool> show_stats_toggle;
};
