/*******************************************
Realization for TARGET class
********************************************
Author: Tsviatkou Anton
RP&E, Belarusian State University 2007
********************************************/
#include "Targ.h"
#include <math.h>
Define_Module(Targ);

Targ::Targ()
{
	newevent =  targmesptr= NULL;	
}

Targ::~Targ()
{// Dispose of dynamically allocated the objects
    cancelAndDelete(newevent);
}
void Targ::initialize(){
	updateInterval = parentModule()->par("updateinterval").doubleValue()/freq;
	speed = par("speed");
	angle = par("angle");
	acceleration = par("acceleration");
	x=par("x");
	y=par("y");
	char display[30];
	sprintf(display, "p=%d,%d;i=status/active", int(x),int(y));
	this->setDisplayString("p=$x,$y;i=status/active");
	s_power=par("power");
	newevent =new cMessage("targevent");
	scheduleAt(updateInterval, newevent);
	counter=0;
}
void Targ::move()
{
	switch (unsigned char(par("moveModel"))){
		case 1://linear motion
			y+=(speed*updateInterval + acceleration*updateInterval*updateInterval/2)*sin(angle*3.14159/180);
			x+=(speed*updateInterval + acceleration*updateInterval*updateInterval/2)*cos(angle*3.14159/180);
			break;
		case 2://linear motion with sinusoid deviations
			y=(1+(sin(2*simTime())/3))*(parentModule()->par("dim_y").doubleValue())/2;
			x+=(speed*updateInterval + acceleration*updateInterval*updateInterval/2)*cos(angle*3.14159/180);
			break;
		case 3://circle motion
			angle += 4*speed/parentModule()->par("dim_x").doubleValue() * updateInterval;
			x = parentModule()->par("dim_x").doubleValue()*(2+cos(angle))/4;
			y = parentModule()->par("dim_x").doubleValue()*(2+sin(angle))/4;
			break;
		default://random motion
			y=normal(y,4,0);
			x+=(speed*updateInterval + acceleration*updateInterval*updateInterval/2)*cos(angle*3.14159/180);
	}
	/*Boundary behavior:
	if object is out of the field it enters from opposite side*/
	if (x>parentModule()->par("dim_x").doubleValue()+20) x=0;
	if (y>parentModule()->par("dim_y").doubleValue()+20) y=0;
	if (x<-20) x=parentModule()->par("dim_x").doubleValue()-x; 
	if(y<-20) y=parentModule()->par("dim_y").doubleValue()-y;
	speed+=acceleration*updateInterval;
	char display[30];
	sprintf(display, "p=%d,%d;i=status/active", int(x),int(y));
	this->setDisplayString(display);
	sprintf(display, "x: %d; y: %d", int(x),int(y));
	displayString().setTagArg("t",0,display);
	par("x")=x;par("y")=y;
	
}
void Targ::handleMessage(cMessage *msg)
{
	this->move();
	if (counter==freq){
		ev << "Generating message from target\n";
		targmesptr = generateMessage();
		ev << targmesptr << endl;
		for (int i=0; i<gate("out")->size()-1; i++) targetsend(targmesptr,i);
		//send message original on last gate
		send(targmesptr, "out", gate("out")->size()-1);
		counter=0;
	}
	counter++;
	//Object display update 'freq' times higher that sensor network tracking update interval
	scheduleAt(simTime()+updateInterval, newevent);
}
Tmes *Targ::generateMessage()
{
    Tmes *msg = new Tmes("Target-to-node",1);//1 means green visualization color of the message
    msg->setS_power(s_power);
	msg->setTarg_x(x);
	msg->setTarg_y(y);
    msg->setTimestamp(simTime());
    return msg;
}
void Targ::targetsend(Tmes *msg, int node)
{
	ev << "TARGET sends message " << msg << " for node " << node << "\n";
	/*Making a copy of the message every time
	This guarantees broadcasting the same message later*/
	cMessage *copy = (cMessage *) msg->dup();
	send(copy, "out", node);   
}
void Targ::finish()
{/*Now we simply remove all connections so we can see
results, which are drawn in Node module*/
for (int i=0; i<gate("out")->size(); i++) gate("out",i)->disconnect();
}

