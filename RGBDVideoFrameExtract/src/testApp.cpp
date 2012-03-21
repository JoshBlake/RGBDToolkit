#include "testApp.h"

//--------------------------------------------------------------
void testApp::setup(){

	moviesLoaded = false;
	ofSystemAlertDialog("Select directory of calibration clips");
	ofFileDialogResult r = ofSystemLoadDialog("Calibration Directory", true);
	if(r.bSuccess){
		movieDirectory = ofDirectory(r.getPath());
		movieDirectory.allowExt("mov");
		movieDirectory.allowExt("mp4");
		movieDirectory.listDir();
		if(movieDirectory.numFiles() > 0){
			for(int i = 0; i < movieDirectory.numFiles(); i++){
				ofVideoPlayer p;
				p.loadMovie( movieDirectory.getPath(i) );
				videoplayers.push_back( p );
				frameExtracted.push_back( false );
			}
			currentMovie = 0;
			moviesLoaded = true;
		}
	}
}

//--------------------------------------------------------------
void testApp::update(){
	if(moviesLoaded){
		videoplayers[currentMovie].setFrame(ofMap(ofGetMouseX(), 0, ofGetWidth(), 0, videoplayers[currentMovie].getTotalNumFrames(), true));
		videoplayers[currentMovie].update();
	}
}

//--------------------------------------------------------------
void testApp::draw(){
	ofBackground(0);
	if(moviesLoaded){
		//calculate a letterboxed preview rectangle based on the source rect compared to the preview rect
		ofRectangle drawRect;
		ofRectangle screenRect = ofRectangle(0,0, ofGetWidth(), ofGetHeight());
		ofRectangle sourceRect = ofRectangle(0,0, videoplayers[currentMovie].getWidth(), videoplayers[currentMovie].getHeight() );
		float screenAspect = screenRect.width/screenRect.height;
		float sourceAspect = sourceRect.width/sourceRect.height;
		if (sourceAspect < screenAspect) {
			drawRect.height = screenRect.height;
			drawRect.width = drawRect.height * sourceAspect;
			drawRect.y = screenRect.y;
			drawRect.x = screenRect.x + screenRect.width / 2 - drawRect.width / 2;
		}
		else{
			drawRect.width = screenRect.width;
			drawRect.height = drawRect.width / sourceAspect;
			drawRect.x = screenRect.x;
			drawRect.y = screenRect.y + screenRect.height / 2 - drawRect.height / 2;
		}
		ofPushStyle();
		if(frameExtracted[currentMovie]){
			ofSetColor(180, 255, 180); //green tint
		}
		videoplayers[currentMovie].draw(drawRect);
		ofPopStyle();
		ofDrawBitmapString("Movie #"+ofToString(currentMovie+1) + " " + movieDirectory.getPath(currentMovie), 10, ofGetHeight()-30);
	}
	else{
		ofSetColor(255);
		ofDrawBitmapString("No Movies Loaded", 10, 10);
	}
}

//--------------------------------------------------------------
void testApp::keyPressed(int key){
	if(key == OF_KEY_LEFT){
		currentMovie--;
		if(currentMovie < 0){
			currentMovie = movieDirectory.numFiles()-1;
		}		
	}
	else if(key == OF_KEY_RIGHT){
		currentMovie++;
		if(currentMovie == movieDirectory.numFiles()){
			currentMovie = 0;
		}
	}
	
	if(key == ' ' && moviesLoaded){
		ofImage frame;
		frame.setFromPixels(videoplayers[currentMovie].getPixelsRef());
		frame.saveImage(ofFilePath::removeExt(movieDirectory.getPath(currentMovie)) + "_calib.png");
		frameExtracted[currentMovie] = true;
	}
}

//--------------------------------------------------------------
void testApp::keyReleased(int key){

}

//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){ 

}