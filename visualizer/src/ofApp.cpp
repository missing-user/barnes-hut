#include "ofApp.h"

void ofApp::initializeParticles() {
  mesh.clear();
  for (const auto &p : particles) {
    mesh.addVertex(p.p);
    auto cur = ofColor::white;
    mesh.addColor(cur);
  }
}
//--------------------------------------------------------------
void ofApp::setup(){
	ofSetVerticalSync(true);
	
  // we're going to load a ton of points into an ofMesh
  mesh.setMode(OF_PRIMITIVE_TRIANGLES);
  ofEnableBlendMode(OF_BLENDMODE_ADD);
  ofSetBackgroundColor(ofColor::black);

	
  particles = make_universe(Distribution::PLUMMER, 1000);
  initializeParticles();

	icoSphere.setMode(OF_PRIMITIVE_TRIANGLES);
	icoSphere.setResolution(0);

  gui.setup();
  gui.add(max_per_node_slider.set("max_per_node", 4, 1, 64));
  gui.add(max_depth_slider.set("max_depth", 32, 1, 48));

  gui.add(timestep_slider.set("timestep", 0.05, 0.001, 0.1));
  //gui.add(brute_force_toggle.set("brute force", false));
  gui.add(theta_slider.set("theta", 1.5, 0.0, 2.5));

  gui.add(show_stats_toggle.set("Show Tree", false));
  gui.add(min_depth_slider.set("min visible level", 4, 0, 20));
}

//--------------------------------------------------------------
void ofApp::update(){
  Tree::maxDepth = max_depth_slider;
  Tree::maxParticles = max_per_node_slider;
	bool brute_force_toggle = false;
  if (brute_force_toggle)
    particles = stepSimulation(particles, timestep_slider);
  else
    particles = stepSimulation(particles, timestep_slider, theta_slider);

 	// Find particle with max v
  auto maxv2 = std::max_element(particles.begin(), particles.end(),
                                [](const Particle &a, const Particle &b) {
                                  return glm::length2(a.v) < glm::length2(b.v);
                                });
  const double maxv_inv = 255.0 / glm::length(maxv2->v);
	// loop through all mesh vertecies and update their positions
  for (std::size_t i = 0; i < mesh.getNumVertices(); i++) {
    mesh.setVertex(i, particles[i].p);
    double len = std::max(50.0, glm::length(particles[i].v) * maxv_inv);
    mesh.setColor(i, ofColor(255 - len / 5, len * 4 / 5, len, 255));
  }
}

void ofApp::drawBoxes(){
	Tree mytree(particles);
  auto boxes = mytree.GetBoundingBoxes();
	auto mdp = mytree.MaxDepthAndParticles();
	//depth_output = ofToString(mdp.first);
	//pcount_output = ofToString(mdp.second);

  ofNoFill();
	for (const auto &b : boxes) {
		if (b.level >= min_depth_slider)
		{
			const auto visualLevel = std::max(0, b.level - min_depth_slider);
			const auto maxLevel = std::max(1, mdp.first - min_depth_slider);
			ofSetColor((255/maxLevel) * visualLevel,
									255 - (255/maxLevel) * visualLevel, 
									0, 128);
			ofDrawBox(b.center, b.dimension.x, b.dimension.y, b.dimension.z);
		}
	}
}

//--------------------------------------------------------------
void ofApp::draw(){
	
  gui.draw();
	ofSetBackgroundColor(ofColor::black);
	ofSetColor(255);
	cam.begin();

	if(show_stats_toggle)
		drawBoxes();
	
	#ifndef TARGET_EMSCRIPTEN
		glPointSize(2);
	#endif
	
  ofFill();
	ofSetColor(ofColor::white, 64);
	// For each particle, draw an icosphere
	icoSphere.setRadius(.1);
	for (const auto &p : particles) {
		icoSphere.setPosition(p.p);
		//icoSphere.setRadius(p.m);
		icoSphere.draw();
	}
	
	cam.end();
	
	// Highlight the nearest vertex
	int n = mesh.getNumVertices();
	float nearestDistance = 0;
	glm::vec2 nearestVertex;
	int nearestIndex = 0;
	glm::vec3 mouse(mouseX, mouseY,0);
	for(int i = 0; i < n; i++) {
		glm::vec3 cur = cam.worldToScreen(mesh.getVertex(i));
		float distance = glm::distance(cur, mouse);
		if(i == 0 || distance < nearestDistance) {
			nearestDistance = distance;
			nearestVertex = cur;
			nearestIndex = i;
		}
	}
	
	ofSetColor(ofColor::gray);
	ofDrawLine(nearestVertex, mouse);
	
  ofNoFill();
	ofSetColor(ofColor::yellow);
	ofSetLineWidth(2);
	ofDrawCircle(nearestVertex, 4);
	ofSetLineWidth(1);
	
	ofDrawBitmapStringHighlight(ofToString(1.0/ofGetLastFrameTime())+" fps", 10,10);
	glm::vec2 offset(10, -10);
	ofDrawBitmapStringHighlight(ofToString(particles.at(nearestIndex).p), mouse + offset);
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
