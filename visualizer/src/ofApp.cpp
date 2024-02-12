#include "ofApp.h"
#include <algorithm>
#include <chrono>

std::vector<Particle> particles;
Particles bodies{64};

void ofApp::initializeParticles() {
  for (size_t i = 0; i <bodies.size(); i++)
  {
    bodies.p[i] = particles[i].p;
    bodies.v[i] = particles[i].v;
    bodies.m[i] = particles[i].m;
  }

  for (size_t i = 0; i <bodies.size(); i++)
  {
    mesh.addVertex(myvec3(bodies.p.x[i], bodies.p.y[i], bodies.p.z[i]));
    auto cur = ofColor::white;
    mesh.addColor(cur);
  }
  std::cout << "Initialized " << bodies.size() << " particles\n";
}

//--------------------------------------------------------------
void ofApp::setup() {
  // ofSetVerticalSync(true);

  // we're going to load a ton of points into an ofMesh
  mesh.setMode(OF_PRIMITIVE_POINTS);
  ofEnableBlendMode(OF_BLENDMODE_ADD);
  ofSetBackgroundColor(ofColor::black);
  glEnable(GL_POINT_SMOOTH); // use circular points instead of square points
  glPointSize(4);            // make the points bigger

  gui.setup();
  gui.add(max_per_node_slider.set("max_per_node", 1, 1, 64));
  gui.add(max_depth_slider.set("max_depth", 64, 1, 96));

  gui.add(timestep_slider.set("timestep", 0.001, 0.0001, 0.01));
  gui.add(brute_force_toggle.set("brute force", true));
  gui.add(theta_slider.set("theta", 1.5, 0.0, 2.5));

  gui.add(num_particles_slider.set("num_particles", 64, 64, 5e4));
  gui.add(mass_slider.set("particle mass", 50, 10.0, 10000.0));
  gui.add(text_output.set("frame time", "text"));

  gui.add(show_stats_toggle.set("Show Tree Stats", false));
  gui.add(min_depth_slider.set("min visible level", 4, 0, 32));
  gui.add(depth_output.set("tree depth", "?"));
  gui.add(pcount_output.set("max particles leafs", "?"));

  particles = make_universe(Distribution::PLUMMER, num_particles_slider);
  initializeParticles();
}
std::vector<DrawableCuboid> drawcuboids;

//--------------------------------------------------------------
void ofApp::update() {
  //Tree::maxDepth = max_depth_slider;
  //Tree::maxParticles = max_per_node_slider;

  set_mass(particles, mass_slider);
  auto begin = std::chrono::steady_clock::now();


    if (brute_force_toggle)
      stepSimulation(bodies, timestep_slider);
    else
      drawcuboids = stepSimulation(bodies, timestep_slider, theta_slider);

  auto end = std::chrono::steady_clock::now();
  auto elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

  text_output = std::to_string(elapsed.count()) + " ms";

  // loop through all mesh vertecies and update their positions
  for (std::size_t i = 0; i < mesh.getNumVertices(); i++) {
    mesh.setVertex(i, myvec3(bodies.p.x[i], bodies.p.y[i], bodies.p.z[i]));
    //double len = std::max(50.0, glm::length(bodies.v[i]) * maxv_inv);
    mesh.setColor(i, ofColor(255 - 255.*i/mesh.getNumVertices(), 255.*i/mesh.getNumVertices(), 255.*i/mesh.getNumVertices(), 255));
  }
}

//--------------------------------------------------------------
void ofApp::draw() {
  gui.draw();
  cam.begin();

  //Tree mytree{bodies};
  //auto boxes = mytree.GetBoundingBoxes();
  ofNoFill();

  if(show_stats_toggle)
  {
    for (const auto &b : drawcuboids) {
      if (b.level >= min_depth_slider)
      {
        const auto visualLevel = std::max(0, b.level - min_depth_slider);
        const auto maxLevel = std::max(1, 16 - min_depth_slider);
        ofSetColor((255/maxLevel) * visualLevel,
                    255 - (255/maxLevel) * visualLevel, 
                    0, 128);
        ofDrawBox(b.center, b.dimension.x, b.dimension.y, b.dimension.z);
      }
    }
    ofSetColor(255);
  }
  mesh.draw();
  cam.end();
  ofDrawBitmapString(
      "   For the \"Big-Bang simulation\" particle distribution example, press R \n\
  For a \"Universe Simulation\" particle distribution example, press T \n \
  To change the particle interaction force to equal the Lennard-Jones molecular Force, press E \n \
  To reorder the array of particles to follow a Morton (Z-order) curve, press O \n \
  To toggle the visualization of the particle order in the stored array, prss L \n \
  To continuously remove a fraction of the kinetic energy of all particles, hold K",
      200, 30);
}

std::vector<Particle> &remove_energy(std::vector<Particle> &particles,
                                     myfloat s = 0.9) {
  for (auto &p : particles) {
    p.v = p.v * s;
  }
  return particles;
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key) {
  if (key == 'r') {
    mesh.clear();
    particles = make_universe(Distribution::BIGBANG, num_particles_slider);
    initializeParticles();
  }
  if (key == 'd') {
    mesh.clear();
    particles = make_universe(Distribution::DEBUG_CUBE, num_particles_slider);
    initializeParticles();
  }
  if (key == 'p') {
    mesh.clear();
    particles = make_universe(Distribution::PLUMMER, num_particles_slider);
    initializeParticles();
  }
  if (key == 't') {
    mesh.clear();
    particles = make_universe(Distribution::UNIVERSE4, num_particles_slider);
    initializeParticles();
  }
  if (key == 'e') {
    mesh.clear();
    particles = make_universe(Distribution::CRYSTALLINE, num_particles_slider);
    initializeParticles();
  }

  if (key == 'o') {
    computeAndOrder(bodies, bounding_box(bodies.p, bodies.size()));
  }

  if (key == 'l') {
    mesh.setMode(mesh.getMode() != OF_PRIMITIVE_LINE_STRIP
                     ? OF_PRIMITIVE_LINE_STRIP
                     : OF_PRIMITIVE_POINTS);
  }

  if (key == 'k') {
    remove_energy(particles);
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
