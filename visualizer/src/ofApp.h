#pragma once

#include "Particle.h"
#include "ofMain.h"
#include "ofxGui.h"
#include <vector>

class ofApp : public ofBaseApp {

private:
  void initializeParticles(int num_particles);

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

  ofEasyCam cam;
  ofMesh mesh;
  ofxPanel gui;

  // Adjustable parameters
  ofxFloatSlider timestep_slider;
  ofxIntSlider max_per_node_slider;
  ofxIntSlider num_particles_slider;
  ofxFloatSlider theta_slider;
  ofxFloatSlider mass_slider;
  ofxToggle brute_force_toggle;

  std::vector<Particle> particles;
  myfloat max_mass;
};
