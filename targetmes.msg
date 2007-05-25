//*******************************************
//Definition for target message called Tmes
//********************************************
//Tmes based on cMessage
//
//Description:
//Target sends this message to nodes of the network
//
//Author: Tsviatkou Anton
//RP&E, Belarusian State University 2007
//********************************************/
message Tmes
{
	fields:
		double s_power;      //Starting intensity of signal
		double timestamp;   //Time when signal was send
		int targ_x;        //x coordinate of the target, to calculate detection error
		int targ_y;  	  //y coordinate of the target, to calculate detection error
}