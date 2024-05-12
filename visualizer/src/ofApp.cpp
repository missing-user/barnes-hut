#include "ofApp.h"
#include <algorithm>
#include <chrono>
#include "Distributions.h"
#include "Barneshut.h"
#include "Cuboid.h"
#include "Order.h"
#include "Particle.h"
#include "Simulation.h"
#include "Forces.h"
Particles particles{};

void ofApp::initializeParticles()
{
  for (size_t i = 0; i < particles.size(); i++)
  {
    mesh.addVertex(myvec3(particles.p.x[i], particles.p.y[i], particles.p.z[i]));
    auto cur = ofColor::white;
    mesh.addColor(cur);
  }
  std::cout << "Initialized " << particles.size() << " particles\n";
}

//--------------------------------------------------------------
void ofApp::setup()
{
  ofSetVerticalSync(true);

  // we're going to load a ton of points into an ofMesh
  mesh.setMode(OF_PRIMITIVE_POINTS);
  ofEnableBlendMode(OF_BLENDMODE_ADD);
  ofSetBackgroundColor(ofColor::black);
  glEnable(GL_POINT_SMOOTH); // use circular points instead of square points
  glPointSize(2);            // make the points bigger

  gui.setup();
  gui.add(max_per_node_slider.set("max_per_node", 1, 1, 64));
  gui.add(max_depth_slider.set("max_depth", 64, 1, 96));

  gui.add(pause_button.set("paused", false));
  gui.add(timestep_slider.set("timestep", 0.001, 0.0001, 0.01));
  gui.add(brute_force_toggle.set("brute force", false));
  gui.add(theta_slider.set("theta", 1.5, 0.0, 2.5));

  gui.add(num_particles_slider.set("num_particles", 256, 10, 5e5));
  gui.add(mass_slider.set("particle mass", 50, 10.0, 1000.0));
  gui.add(text_output.set("frame time", "text"));

  gui.add(show_stats_toggle.set("Show Tree Stats", false));
  gui.add(min_depth_slider.set("min visible level", 4, 0, 32));
  gui.add(depth_output.set("tree depth", "?"));
  gui.add(pcount_output.set("max particles leafs", "?"));

  particles = make_universe(Distribution::PLUMMER, num_particles_slider);
  initializeParticles();
}

//--------------------------------------------------------------
void ofApp::update()
{
  // Tree::maxDepth = max_depth_slider;
  // Tree::maxParticles = max_per_node_slider;
  if (!pause_button)
  {
    set_mass(particles, mass_slider);
    auto begin = std::chrono::steady_clock::now();

    if (brute_force_toggle)
      stepSimulation(particles, timestep_slider);
    else
      stepSimulation(particles, timestep_slider, theta_slider);

    auto end = std::chrono::steady_clock::now();
    auto elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

    text_output = std::to_string(elapsed.count()) + " ms";
  }
}

//--------------------------------------------------------------
void ofApp::draw()
{
  gui.draw();
  cam.begin();
  ofNoFill();

  // Find maximumum velocity
  double vmax = 0;
  for (std::size_t i = 0; i < mesh.getNumVertices(); i++)
  {
    vmax = std::max(glm::length(myvec3(particles.v.x[i], particles.v.y[i], particles.v.z[i])), vmax);
  }

  // loop through all mesh vertecies and update their positions
  for (std::size_t i = 0; i < mesh.getNumVertices(); i++)
  {
    mesh.setVertex(i, myvec3(particles.p.x[i], particles.p.y[i], particles.p.z[i]));
    double particle_v = glm::length(myvec3(particles.v.x[i], particles.v.y[i], particles.v.z[i]));
    mesh.setColor(i, ofColor(255 - 255. * particle_v / vmax, 255. * particle_v / vmax, 255. * particle_v / vmax, 255));
  }

  if (show_stats_toggle)
  {
    std::cout << particles.size() << std::endl;
    std::cout << particles.v.y[5] << std::endl;
    myvec3 evalpos = {100, 0, 0};
    auto info = bh_superstep_debug(evalpos, particles, theta_slider * theta_slider);
    depth_output = std::to_string(info.depth);
    pcount_output = std::to_string(info.max_particles_in_leaf);

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

Particles &remove_energy(Particles &particles,
                         myfloat s = 0.9)
{
  for (auto &x : particles.v.x)
  {
    x *= s;
  }
  for (auto &y : particles.v.y)
  {
    y *= s;
  }
  for (auto &z : particles.v.z)
  {
    z *= s;
  }
  return particles;
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key)
{
  if (key == 'r')
  {
    mesh.clear();
    particles = make_universe(Distribution::BIGBANG, num_particles_slider);
    initializeParticles();
  }
  if (key == 'd')
  {
    mesh.clear();
    particles = make_universe(Distribution::DEBUG_CUBE, num_particles_slider);
    initializeParticles();
  }
  if (key == 'p')
  {
    mesh.clear();
    particles = make_universe(Distribution::PLUMMER, num_particles_slider);
    initializeParticles();
  }
  if (key == 't')
  {
    mesh.clear();
    particles = make_universe(Distribution::UNIVERSE4, num_particles_slider);
    initializeParticles();
  }
  if (key == 'e')
  {
    mesh.clear();
    particles = make_universe(Distribution::CRYSTALLINE, num_particles_slider);
    initializeParticles();
  }

  if (key == 'o')
  {
    computeAndOrder(particles, bounding_box(particles.p));
  }

  if (key == 'l')
  {
    mesh.setMode(mesh.getMode() != OF_PRIMITIVE_LINE_STRIP
                     ? OF_PRIMITIVE_LINE_STRIP
                     : OF_PRIMITIVE_POINTS);
  }

  if (key == 'k')
  {
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
