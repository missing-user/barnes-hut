#include "ofApp.h"
#include "Distributions.h"
#include "Simulation.h"
#include <chrono>

void ofApp::initializeParticles() {
  max_mass = 0;

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
  // ofSetVerticalSync(true);

  // we're going to load a ton of points into an ofMesh
  mesh.setMode(OF_PRIMITIVE_POINTS);

  glEnable(GL_POINT_SMOOTH); // use circular points instead of square points
  glPointSize(2);            // make the points bigger

  calcDepthButton.addListener(this, &ofApp::calcDepthButtonPressed);

  gui.setup();
  gui.add(max_per_node_slider.set("max_per_node", 1, 1, 128));
  gui.add(max_depth_slider.set("max_depth", 32, 1, 128));

  gui.add(timestep_slider.set("timestep", 0.001, 0.001, 0.1));
  gui.add(brute_force_toggle.set("brute force", false));
  gui.add(theta_slider.set("theta", 1.5, 0.0, 2.5));

  gui.add(num_particles_slider.set("num_particles", 1000, 100, 1e5));
  gui.add(mass_slider.set("particle mass", 50, 10.0, 10000.0));
  gui.add(text_output.set("frame time", "text"));

  gui.add(calcDepthButton.setup("Compute actual tree depth"));
  gui.add(depth_output.set("tree depth", "?"));
  gui.add(pcount_output.set("max particles leafs", "?"));

  particles = make_universe(Distribution::BIGBANG, num_particles_slider);
  initializeParticles();
}

//--------------------------------------------------------------
void ofApp::update() {
  Tree::maxDepth = max_depth_slider;
  Tree::maxParticles = max_per_node_slider;

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
    particles = make_universe(Distribution::BIGBANG, num_particles_slider);
    initializeParticles();
  }
  if (key == 't') {
    mesh.clear();
    particles = make_universe(Distribution::UNIVERSE4, num_particles_slider);
    initializeParticles();
  }
  if (key == 'e') {
    mesh.clear();
    particles = make_universe(Distribution::COLLISION, num_particles_slider);
    initializeParticles();
  }
}

void ofApp::calcDepthButtonPressed() {
  Tree mytree{particles};
  auto minmax = mytree.MaxDepthAndParticles();
  depth_output = std::to_string(minmax.first);
  pcount_output = std::to_string(minmax.second);
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
