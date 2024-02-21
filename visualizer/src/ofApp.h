#pragma once

#include "ofMain.h"
#include "ofxGui.h"
#include <string>
#include <vector>

class ofApp : public ofBaseApp
{

private:
  void initializeParticles();

  ofEasyCam cam;
  ofMesh mesh;
  ofxPanel gui;

  // Adjustable parameters
  ofParameter<double> timestep_slider;
  ofParameter<int> max_depth_slider;
  ofParameter<int> max_per_node_slider;
  ofParameter<int> num_particles_slider;
  ofParameter<double> theta_slider;
  ofParameter<double> mass_slider;
  ofParameter<bool> brute_force_toggle;

  // text output
  ofParameter<std::string> text_output;

  // stats
  ofParameter<std::string> depth_output;
  ofParameter<std::string> pcount_output;
  ofParameter<int> min_depth_slider;
  ofParameter<bool> show_stats_toggle;

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
  void HUDGuide();
};
