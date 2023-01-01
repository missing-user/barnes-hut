#include "ofApp.h"
#include "Distributions.h"
#include "Simulation.h"
#include <chrono>

void ofApp::initializeParticles(int num_particles = 1000) {
  max_mass = 0;
  particles = bigbang(num_particles); // stable_orbit(); // universe1();
  for (const auto &p : particles) {
    if (p.m > max_mass) {
      max_mass = p.m;
    }
  }

  for (const auto &p : particles) {
    mesh.addVertex(p.p);
    auto cur = ofColor::white;
    cur.b = (p.m / max_mass) * 255;
    mesh.addColor(cur);
  }
}

//--------------------------------------------------------------
void ofApp::setup() {

  ofSetVerticalSync(true);

  // we're going to load a ton of points into an ofMesh
  mesh.setMode(OF_PRIMITIVE_POINTS);

  initializeParticles();

  // ofEnableDepthTest();
  glEnable(GL_POINT_SMOOTH); // use circular points instead of square points
  glPointSize(2);            // make the points bigger

  gui.setup();
  gui.add(timestep_slider.setup("timestep", 0.001, 0.001, 0.1));
  gui.add(max_per_node_slider.setup("max_per_node", 32, 1, 128));
  gui.add(num_particles_slider.setup("num_particles", 1000, 100, 10000));
  gui.add(theta_slider.setup("theta", 1.5, 0.0, 2.5));
  gui.add(mass_slider.setup("particle mass", 50, 1.0, 1000.0));
  gui.add(brute_force_toggle.setup("brute force", false));
}

//--------------------------------------------------------------
void ofApp::update() {

  auto begin = std::chrono::steady_clock::now();
  Tree::setMaxPointsPerNode(max_per_node_slider);
  set_mass(particles, mass_slider);

  if (!brute_force_toggle)
    particles = stepSimulation(particles, timestep_slider, theta_slider);
  else
    particles = stepSimulation(particles, timestep_slider);

  auto end = std::chrono::steady_clock::now();
  auto elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
  std::cout << "Simulation frame time: " << elapsed.count() << " ms"
            << std::endl;

  // loop through all mesh vertecies and update their positions
  for (int i = 0; i < mesh.getNumVertices(); i++) {
    mesh.setVertex(i, particles[i].p);
  }
}

//--------------------------------------------------------------
void ofApp::draw() {
  ofBackgroundGradient(ofColor::gray, ofColor::black, OF_GRADIENT_CIRCULAR);
  gui.draw();

  cam.begin();
  ofRotateYDeg(90);
  // ofTranslate(-img.getWidth() / 2, -img.getHeight() / 2);
  mesh.draw();
  cam.end();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key) {
  if (key == 'r') {
    mesh.clear();
    initializeParticles(num_particles_slider);
  }
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y) {}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button) {}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button) {}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button) {}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y) {}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y) {}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h) {}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg) {}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo) {}
