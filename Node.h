/*******************************************
Definition for NODE class
********************************************
Author: Tsviatkou Anton
RP&E, Belarusian State University 2007
********************************************/
#include <stdio.h>
#include <string.h>
#include <omnetpp.h>
#include "mes_m.h"
#include "targetmes_m.h"
#include "node-master_m.h"
/*DEFINITION
1.Information structure
Used in NODE messages when alerting the others
2.Target coordinates structure
3.Structure to keep the track
*/
const int polypower=3;
const int appr_points=10;
struct information{
double intensity; //node intensity reading
double x,y;		  //Coordinates x,y
};
struct targ_coords{
double x;
double y;
};
struct track{
double intensity;
double x;
double y;
double targ_x;
double targ_y;
double err_x;
double err_y;
};
std::map<double,track> m_track; //data container with calculated coordinates
typedef std::multimap<double,information>::const_iterator iter;
typedef std::map<double,track>::iterator trackiter;
class Node : public cSimpleModule
{  protected:
	cOutVector trajectoryX;
	cOutVector trajectoryY;
	cOutVector trajectoryapprY;
    virtual Nodemes *generatePacket(double timestamp, information targinfo);//generates message packet
    virtual void initialize();
    virtual void handleMessage(cMessage *msg);//called when message arrives
    virtual void record(double timestamp, information nodedata);//records to 'data' multimap container
	virtual targ_coords calculate(double timestamp);//realizes triangulation algorithm
	virtual void send_packet(cMessage *msg, int sendgate);//used to send a copy of the same message for all neighbours
	virtual void finish();//recording statistics & display estimation errors
	virtual void Node::sendToMaster(targ_coords coords,double timestamp, double intensity);
	double Node::calculateApproximationError(double time0, double time1, double max_dy,bool flag);
	bool Node::setLinearInterpolationCoefficients(double x[appr_points], double y[appr_points]);
public:
	Node();
	virtual ~Node();
	void SetCalculatedTime(double time);
private:
    cMessage *packetbroadcast;//simulator event message which starts node broadcast 
	cMessage *packet;//node sends it to alert the neighbours
	double threshold;//intensity threshold, if measurment>threshold node sees the target
	double time;//to keep timestamp of the last calculation(helps to prevent recalculation)
	long numSent;//num times of transmission for other nodes
    long numTReceived;//num messages from target
    long numNReceived;//num messages from nodes
	long numCalc;//num times make the calculation of target position
	std::multimap<double,information> data; //NODE data container with adjasent nodes measurments and coordinates
};