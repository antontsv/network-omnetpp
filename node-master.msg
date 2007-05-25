//*******************************************
//Definition for target message called ForMastermes
//********************************************
//ForMastermes based on cMessage
//
//Description:
//After target position calculation leader node 
//will send this message to masternode through the network
//
//Author: Tsviatkou Anton
//RP&E, Belarusian State University 2007
//********************************************/
message ForMastermes
{
	fields:
		int nodeID;                //Leader node ID
		double timestamp;         //information timestamp
		double targ_x;           //calculated x coordinate of the target
		double targ_y;  	    //calculated y coordinate of the target
		double distance;       //used for shortest path routing simulation
		double intensity;     //receved target's signal intensity
}