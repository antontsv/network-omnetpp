/*******************************************
Definition for TARGET class
********************************************
Author: Tsviatkou Anton
RP&E, Belarusian State University 2007
********************************************/
#include "targetmes_m.h"
const unsigned char freq=10;
class Targ : public cSimpleModule
{protected:
	virtual void initialize();
    virtual void move();
	virtual Tmes *generateMessage();
    virtual void targetsend(Tmes *msg, int node);
    virtual void handleMessage(cMessage *msg);
	virtual void finish();
public:
	Targ();
	virtual ~Targ();	
private:
	cMessage *newevent;
	Tmes *targmesptr;
	double x;
	double y;
	double speed;
	double angle; //range [0,360)
	double acceleration;
	double updateInterval;
	int s_power; //target signal power
	unsigned char counter;//used for target position on-screen update
};