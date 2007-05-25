//*******************************************
//Definition for node message called Nodemes
//********************************************
//Nodemes based on cMessage
//
//Description:
//Node sends this message to other nodes when receive signal from target
//i.e. Nodemes is sent upon target detection to help other nodes to detect
//target coordinates
//
//Author: Tsviatkou Anton
//RP&E, Belarusian State University 2007
//********************************************/
message Nodemes
{
	fields:
		double targintensity; //Intensity of received signal
		double timestamp;     //Time when signal was send by target
		double timedelay;     //arrival delay time
		int node_x;           //x coordinate of the node, that sends Nodemes
		int node_y;           //y coordinate of the node, that sends Nodemes
}