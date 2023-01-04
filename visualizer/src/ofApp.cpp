#include "ofApp.h"
#include "Distributions.h"
#include "Simulation.h"
#include <chrono>

void ofApp::initializeParticles(int num_particles = 1000, bool type = false) {
  max_mass = 0;
  if (type)
    particles = universe4(num_particles);
  else
    particles = bigbang(num_particles);

  for (const auto &p : particles) {
    if (p.m > max_mass) {
      max_mass = p.m;
    }
  }

  for (const auto &p : particles) {
    mesh.addVertex(p.p);
    auto cur = ofColor::white;
    cur = ofColor(ofRandom(128, 255), ofRandom(128, 255), ofRandom(128, 255));

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

  glEnable(GL_POINT_SMOOTH); // use circular points instead of square points
  glPointSize(3);            // make the points bigger

  gui.setup();
  gui.add(timestep_slider.set("timestep", 0.001, 0.001, 0.1));
  gui.add(max_per_node_slider.set("max_per_node", 32, 1, 128));
  gui.add(num_particles_slider.set("num_particles", 1000, 100, 10000));
  gui.add(theta_slider.set("theta", 1.5, 0.0, 2.5));
  gui.add(mass_slider.set("particle mass", 50, 10.0, 10000.0));
  gui.add(brute_force_toggle.set("brute force", false));
  gui.add(text_output.set("frame time", "text"));
}

//--------------------------------------------------------------
void ofApp::update() {
  Tree::setMaxDepth(max_per_node_slider);

  set_mass(particles, mass_slider);
  auto begin = std::chrono::steady_clock::now();

  if (!brute_force_toggle)
    particles = stepSimulation(particles, timestep_slider, theta_slider);
  else
    particles = stepSimulation(particles, timestep_slider);

  auto end = std::chrono::steady_clock::now();
  auto elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

  text_output = std::to_string(elapsed.count()) + " ms";

  // loop through all mesh vertecies and update their positions
  for (int i = 0; i < mesh.getNumVertices(); i++) {
    mesh.setVertex(i, particles[i].p);
    double len = std::max(10.0, std::min(glm::length(particles[i].v), 255.0));
    mesh.setColor(i, ofColor(255, len, len));
  }
}

//--------------------------------------------------------------
void ofApp::draw() {
  ofBackgroundGradient(ofColor::black, ofColor::black, OF_GRADIENT_CIRCULAR);
  gui.draw();

  cam.begin();
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
  if (key == 't') {
    mesh.clear();
    initializeParticles(num_particles_slider, true);
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
