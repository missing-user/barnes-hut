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
	mesh.load("lofi-bunny.ply");

	
  // we're going to load a ton of points into an ofMesh
  //mesh.setMode(OF_PRIMITIVE_POINTS);
  ofEnableBlendMode(OF_BLENDMODE_ADD);
  ofSetBackgroundColor(ofColor::black);

	
  particles = make_universe(Distribution::BIGBANG, 1000);
  initializeParticles();
}

//--------------------------------------------------------------
void ofApp::update(){

}

void ofApp::drawBoxes(){
	int min_depth_slider = 3;

	Tree mytree(particles);
  auto boxes = mytree.GetBoundingBoxes();
  ofNoFill();
    auto mdp = mytree.MaxDepthAndParticles();
    //depth_output = ofToString(mdp.first);
    //pcount_output = ofToString(mdp.second);

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
    ofSetColor(255);
}

//--------------------------------------------------------------
void ofApp::draw(){
	//ofBackgroundGradient(ofColor(64), ofColor(0));
	
	ofSetColor(255);
	cam.begin();

	drawBoxes();
	
	ofSetColor(ofColor::gray);
	mesh.drawWireframe();
	
	#ifndef TARGET_EMSCRIPTEN
		glPointSize(2);
	#endif
	ofSetColor(ofColor::white);
	mesh.drawVertices();
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
	
	ofSetColor(ofColor::yellow);
	ofSetLineWidth(2);
	ofDrawCircle(nearestVertex, 4);
	ofSetLineWidth(1);
	
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
