/*******************************************
Realization for NODE class
********************************************
Author: Tsviatkou Anton
RP&E, Belarusian State University 2007
********************************************/
#include "Node.h"
#include <math.h>

int max_elements=10;//NODE data base volume for temporary readings
unsigned char masterNodeID;
int masterNodeTrasmissions=0;
Define_Module(Node);//Define classes with OMNeT++ simulator
double C[4];

/*Parameters for error estimation at the end*/
double tx;
double ty;
double speed;
double angle;
double acceleration;
int moveModel;

Node::Node()
{
packetbroadcast = packet  =  NULL;
numSent = numTReceived = numNReceived = numCalc =0;
time=0;
}

Node::~Node()
{// Dispose of dynamically allocated the objects
cancelAndDelete(packetbroadcast);
}
void Node::initialize()
{
	if (index()==0){/*Parameters initiazation. They will be used for error estimation at the end*/
	tx=parentModule()->submodule("target")->par("x").doubleValue();
	ty=parentModule()->submodule("target")->par("y").doubleValue();
	speed=parentModule()->submodule("target")->par("speed").doubleValue();
	angle=parentModule()->submodule("target")->par("angle").doubleValue();
	acceleration=parentModule()->submodule("target")->par("acceleration").doubleValue();
	moveModel=parentModule()->submodule("target")->par("moveModel");
	ev<<"INITIAL:\nx="<<tx<<" y="<<ty<<" speed="<<speed<<" angle="<<angle<<"\n";
	ev<<"acceleration="<<acceleration<<" moveModel="<<moveModel<<"\n";
	}
	trajectoryX.setName("x");
	trajectoryY.setName("y");
	trajectoryapprY.setName("appry");
	packetbroadcast=new cMessage("Broadcast");
	threshold=par("threshold");
	WATCH(numSent);
	WATCH(numTReceived);
	WATCH(numNReceived);
	WATCH(numCalc);
	int num=parentModule()->par("numNodes");//get total number of nodes
	int numl=int(sqrt(double(num)));
	int dim_x=parentModule()->par("dim_x"),dim_y=parentModule()->par("dim_y");
	int dx=dim_x/(numl-1),dy=dim_y/numl;
	/*Now we define node coordinates assuming matrix architecture*/
	int y=50+dy*(index()/numl);
	int x=50+dx*(index()%numl);
	ev<<"Node["<<index()<<"]: x="<<x<<" y="<<y<<"\n";//Display coords in the window
	this->par("x")=x;
	this->par("y")=y;
	char display[30];
	if (index()==int(num/2+0.5)){//Define masternode appearance for visualisazion layer
		sprintf(display, "p=%d,%d;i=status/up", x,y);
		masterNodeID=index();}
		else sprintf(display, "p=%d,%d;i=status/off", x,y);//Define node appearance for visualisazion layer
	this->setDisplayString(display);
	if(index()<num-1) return;//if last one rerurn
	/*Now we define how to connect our nodes*/
	for(int i=0;i<num;i++){
		cModule* module=parentModule()->submodule("node",i);//gets pointer to other nodes
		for(int j=0;j<module->gate("out")->size();j++)
		{
			cModule* submodule=module->gate("out",j)->destinationGate()->ownerModule();
			dim_x=submodule->par("x");
			dim_y=submodule->par("y");
			x=module->par("x");
			y=module->par("y");
			if(int(module->par("transmit"))==-1) numl=max(dx,dy); 
				else numl=module->par("transmit");//if transmit power is set in omnetpp.ini
			if((x-dim_x)*(x-dim_x)+(y-dim_y)*(y-dim_y)>numl*numl) module->gate("out",j)->disconnect();
			//in fact all nodes are connected, we disconnect some that do not meet distance-power criteria
		};
	};
}
void Node::record(double timestamp, information nodedata)
{	
	//Now we erase some temporary data to get some storage
	if (time<timestamp) data.erase(data.begin(),data.upper_bound(time));//Erase if data is too old
	if (data.size()>=max_elements) data.erase(--data.end());
	data.insert(std::make_pair(timestamp,nodedata));
	ev<<"#Keys:"<<data.count(timestamp)<<"("<<data.size()<<")\n";
	if (data.count(timestamp)>=3) {//Now we are ready to use triangulation ALG
		targ_coords coords=calculate(timestamp);
	}
}
targ_coords Node::calculate(double timestamp)
{ /*Calculate() gives target coordinates
	by solving a system with 3 eguations.
	It also gives initial target radiation power*/
	std::pair<iter,iter> b=data.equal_range(timestamp);
	//trackiter TI=m_track.find(timestamp);
	information x[3];
	targ_coords coords={-1,-1};	
	if(time==timestamp) {data.erase(data.begin(),data.upper_bound(timestamp));return coords;}//Recalculation attempt check
	int j=0;
	for(iter i=b.first; i!=b.second;++i) {x[j]=(*i).second;j++;}//extract data for ALG
	//ALG beginning
	double a12=2*(x[0].x-x[1].x),a13=2*(x[0].x-x[2].x),b12=2*(x[0].y-x[1].y),b13=2*(x[0].y-x[2].y);
	double g21=x[1].x*x[1].x+x[1].y*x[1].y-x[0].x*x[0].x-x[0].y*x[0].y,g31=x[2].x*x[2].x+x[2].y*x[2].y-x[0].x*x[0].x-x[0].y*x[0].y;
	double d13=(x[0].intensity-x[2].intensity)/(x[0].intensity*x[2].intensity),d12=(x[0].intensity-x[1].intensity)/(x[0].intensity*x[1].intensity);
	double det=a12*b13-a13*b12;
	if (det==0) {ev<<"CALCULATION forbitten: DET=0\n";data.erase(--data.end());return coords;}//DET check, if true erasing one component
	time=timestamp;//update time the node had a chanse to start ALG
	data.erase(data.begin(),data.upper_bound(timestamp));//erase data from DB, all the data in the buffer
	/*Check if meets our main criteria: calculations should be done by target closest sensor*/
	if(!(x[0].intensity>=x[1].intensity&&x[0].intensity>=x[2].intensity)) {ev<<"CALCULATION forbitten: node intensity is lower than others\n";return coords;}
	double k1=d12*b13-d13*b12,k2=g21*b13-g31*b12+x[1].x*det,l1=d13*a12-d12*a13,l2=g31*a12-g21*a13+x[1].y*det;
	double ksi=det*det/x[1].intensity;
	double A=k1*k1+l1*l1,B=2*k1*k2+2*l1*l2+ksi,C=k2*k2+l2*l2;
	double p0;
	if (A>1e-9){p0=(B-sqrt(abs(B*B-4*A*C)))/(2*A);}
		else p0=-C/B;
	double def_x=((p0*d12-g21)*b13-(p0*d13-g31)*b12)/det;
	double def_y=((p0*d13-g31)*a12-(p0*d12-g21)*a13)/det;
	coords.x=def_x;coords.y=def_y;
	
	//Error estimation
	//if (TI->second.intensity<x[0].intensity){
	//TI->second.intensity=x[0].intensity;
	//TI->second.x=def_x;TI->second.y=def_y;TI->second.err_x=abs(TI->second.targ_x-def_x);TI->second.err_y=abs(TI->second.targ_y-def_y);
	ev<<"Node "<<index()<<":\n";
	ev<<"CALCULATION Info:\n"<<"Def Info: x="<<def_x<<" y="<<def_y<<"\n";
	//ev<<"Err Info: x="<<TI->second.err_x<<" y="<<TI->second.err_y<<"\n";
	//}
	numCalc++;
	sendToMaster(coords, time, x[0].intensity);
	return coords;
}
void Node::handleMessage(cMessage *msg)
{
	if(std::string(msg->name())=="PacketCopy"){/*here we receive message from other nodes*/
		ev<<"Node receives message from node...\n";
		Nodemes *nodemsg = check_and_cast<Nodemes *>(msg);
		information i;
		i.intensity=nodemsg->getTargintensity();
		i.x=nodemsg->getNode_x();
		i.y=nodemsg->getNode_y();
		record(nodemsg->getTimestamp(),i);//saving arrived data to node DB
		delete msg;
		numNReceived++;
		return;
	};
	if(std::string(msg->name())=="MessageForMaster"){/*here we receive message for masternode from other nodes*/
		ForMastermes *nodemsg = check_and_cast<ForMastermes *>(msg);
		if(index()==masterNodeID){
			targ_coords i;
			i.x=nodemsg->getTarg_x();
			i.y=nodemsg->getTarg_y();
			trackiter TI=m_track.find(nodemsg->getTimestamp());
			if (TI->second.intensity<nodemsg->getIntensity()){
					TI->second.intensity=nodemsg->getIntensity();
					TI->second.x=i.x;TI->second.y=i.y;
					TI->second.err_x=100*abs(TI->second.targ_x-i.x)/TI->second.targ_x;
					TI->second.err_y=100*abs(TI->second.targ_y-i.y)/TI->second.targ_y;
			}
		}
		else { 
			int x=int(parentModule()->submodule("node",masterNodeID)->par("x"))-int(this->par("x")),
				y=int(parentModule()->submodule("node",masterNodeID)->par("y"))-int(this->par("y"));
			if (x*x+y*y<nodemsg->getDistance()){
				nodemsg->setDistance(x*x+y*y);
				for (int i=0; i<gate("out")->size(); i++) 
					if(gate("out",i)->isConnected())send((cMessage *) nodemsg->dup(), "out", i);
				numSent++;
			}
		}
		delete msg;
		numNReceived++;
		return;
	};
	if(msg==packetbroadcast){
		bubble("See target!");
		/*node broadcasts a message to alert it's neighbours*/
		ev<<"Broadcasting...\n";
		for (int i=0; i<gate("out")->size(); i++) if(gate("out",i)->isConnected()) send_packet(packet,i);
		delete packet; 
		numSent++;
		return;
	}
	Tmes *targmsg = check_and_cast<Tmes *>(msg);
	double x=par("x");
	double y=par("y");
	x-=targmsg->getTarg_x();
	y-=targmsg->getTarg_y();
	//Check 'above the head' situation
	if((abs(x)<1e-6)&&(abs(y)<1e-6)){
		track AboveHead={0,0,0,targmsg->getTarg_x(),targmsg->getTarg_y(),0,0};
		m_track.insert(std::make_pair(targmsg->getTimestamp(),AboveHead));
		targ_coords coords={targmsg->getTarg_x(),targmsg->getTarg_y()};
		bubble("ABOVE THE HEAD!");
		sendToMaster(coords, targmsg->getTimestamp(), 1000);
		numCalc++;
		int num=parentModule()->par("numNodes");
		cModule *nodeModule;
		Node *currentNode;
		for(int i=0;i<num;i++){
			nodeModule=parentModule()->submodule("node",i);
			currentNode=check_and_cast<Node *>(nodeModule);
			currentNode->SetCalculatedTime(targmsg->getTimestamp());
		}
		delete msg;
		return;
	}
	double intensity=targmsg->getS_power()/(4*3.14159*(x*x+y*y));//MODEL: estimating intensity we should receive
	intensity*=1-0.05*normal(0,1,0);//Adding error to intensity
	ev << "ARRIVED from TARGET: intensity=" <<intensity<<"; was send at "<<targmsg->getTimestamp()<<"\n";
	ev<<"Target real coordinates X="<<targmsg->getTarg_x()<<" Y="<<targmsg->getTarg_y()<<"\n";
	/*If threshold is passed node begins tracking*/
	if (intensity>=threshold){
		scheduleAt(simTime()+parentModule()->par("updateinterval").doubleValue()*(1-1/(x*x+y*y)), packetbroadcast);
		numTReceived++;
		information f;
		f.intensity=intensity;
		f.x=par("x");f.y=par("y");
		track newtrack={0,0,0,targmsg->getTarg_x(),targmsg->getTarg_y(),0,0};
		m_track.insert(std::make_pair(targmsg->getTimestamp(),newtrack));
		record(targmsg->getTimestamp(),f);
		packet = generatePacket(targmsg->getTimestamp(),f);
	}
	delete msg;
}

Nodemes *Node::generatePacket(double timestamp, information targinfo)
{/*Node creates packet to alert about target
 Packet contain important info: node coordinates, node measurements*/
	int src = index();
    char msgname[20];
    sprintf(msgname, "Node[%d]-broadcast", src);
	Nodemes *msg = new Nodemes(msgname);
    msg->setTimestamp(timestamp);
    msg->setTargintensity(targinfo.intensity);
	msg->setNode_x(targinfo.x);
	msg->setNode_y(targinfo.y);
	packet = msg;
    return msg;
}
void Node::send_packet(cMessage *msg, int sendgate)
{//Making a copy of the message every time
//This guarantees broadcasting the same message later
	cMessage *copy = (cMessage *) msg->dup();
	copy->setName("PacketCopy");
	send(copy, "out", sendgate);
}
void Node::SetCalculatedTime(double time){
	this->time=time;
}
void Node::sendToMaster(targ_coords coords, double timestamp, double intensity){
	if (index()!=masterNodeID){
		int x=int(parentModule()->submodule("node",masterNodeID)->par("x"))-int(this->par("x")),
			y=int(parentModule()->submodule("node",masterNodeID)->par("y"))-int(this->par("y"));	
		ForMastermes *msg = new ForMastermes("MessageForMaster",2);//2 means blue visualization color of the message
		msg->setNodeID(index());
		msg->setTimestamp(timestamp);
		msg->setTarg_x(coords.x);
		msg->setTarg_y(coords.y);
		msg->setDistance(x*x+y*y);
		msg->setIntensity(intensity);
		for (int i=0; i<gate("out")->size(); i++) 
			if(gate("out",i)->isConnected())send((cMessage *) msg->dup(), "out", i);
		delete msg;
		numSent++;
	}
	else {
		trackiter TI=m_track.find(timestamp);	
		TI->second.intensity=intensity;
		TI->second.x=coords.x;TI->second.y=coords.y;
		TI->second.err_x=100*abs(TI->second.targ_x-coords.x)/TI->second.targ_x;
		TI->second.err_y=100*abs(TI->second.targ_y-coords.y)/TI->second.targ_y;
	}
return;
}
double power(double x, unsigned char degree){
double result=1;
if (degree==0) return result;
for (unsigned char i=1;i<=degree;i++) result*=x;
return result;
}
int setApproximationCoefficients(double x[appr_points], double y[appr_points]){
	double B[polypower+1],temp;
	double A[polypower+1][polypower+1];
	for(int i=0;i<=polypower;i++){
		B[i]=0;
		for(int j=0;j<=polypower;j++) A[i][j]=0;
	}
	for(int i=0;i<=polypower;i++)
		for(int j=0;j<appr_points;j++)B[i]+=y[j]*power(x[j],i);
	
	/*double POWERX[2*polypower];
	for(int i=0;i<2*polypower;i++){ POWERX[i]=0;
		for(int j=0;j<appr_points;j++) POWERX[i]+=power(x[j],i+1);
	}
	for(int i=0;i<appr_points;i++)
		for(int j=0;j<appr_points;j++) A[i][j]=POWERX[abs(i+j-1)];
	A[0][0]=appr_points;*/
	
	for(int i=0;i<=polypower;i++)
		for(int j=0;j<=polypower;j++)
			for(int k=0;k<appr_points;k++)A[i][j]+=power(x[k],i+j);
	
	//Solving linear system via Gauss method
	unsigned char IOR[polypower+1],l,M,p;//IOR keeps indices of strings in our system
	for(int i=0;i<=polypower;i++) IOR[i]=i;
	for(int k=0;k<polypower;k++){
		temp=0;//now we use it to find max element
		for(int j=k;j<=polypower;j++){
			l=IOR[j];
			if(abs(A[l][k])<temp) continue;
			M=l;p=j;temp=abs(A[l][k]);
		}
		if(p!=k){IOR[p]=IOR[k];IOR[k]=M;}
		temp=A[M][k];//now temp is selected element
		if (temp==0) return 1;//Error occured, returning exit code #1
		for(int j=k;j<=polypower;j++) A[M][j]/=temp;
		B[M]/=temp;
		for(int i=k+1;i<=polypower;i++){
			l=IOR[i];
			for(int j=k+1;j<=polypower;j++) A[l][j]-=A[l][k]*A[M][j];
			B[l]-=A[l][k]*B[M];
			A[l][k]=0;
		}	
	}
	l=IOR[polypower];
	if (A[l][polypower]==0) return 1;//Error occured, returning exit code #1
	B[l]/=A[l][polypower];
	C[polypower]=B[l];
	for(int k=polypower-1;k>=0;k--){
		l=IOR[k];
		C[k]=B[l];
		for(int j=k+1;j<=polypower;j++)	C[k]-=C[j]*A[l][j];
	}
temp=0;
for(int j=0;j<appr_points;j++)temp+=power(y[j]-C[0]-C[1]*x[j]-C[2]*power(x[j],2)-C[3]*power(x[j],3),2);
temp/=appr_points-polypower-1;
//ev<<"DISP: "<<temp<<"  STDEV: "<<sqrt(temp)<<"\n";
	return 0;//Returning code #0, i.e. OK
}
bool Node::setLinearInterpolationCoefficients(double x[appr_points], double y[appr_points]){
if (x[0]==x[1]) return false;
if (abs(1/(x[0]-x[1]))>parentModule()->par("dim_x").doubleValue()) return false;
C[1]=(y[0]-y[1])/(x[0]-x[1]);
if(x[0]!=0)C[0]=y[0]-C[1]*x[0];
	else C[0]=y[1]-C[1]*x[1];
C[2]=0;
C[3]=0;
return true;
}
double avg_err=0,avg_msr=0;
double Node::calculateApproximationError(double time0, double time1, double max_dy,bool flag){
	double updateInterval = 0.0001;//we can use: parentModule()->par("updateinterval").doubleValue()/100;
	double appry;
		for (double i=time0+updateInterval;i<=time1;i+=updateInterval){
			switch (moveModel){
			case 1://linear motion
				ty+=(speed*updateInterval + acceleration*updateInterval*updateInterval/2)*sin(angle*3.14159/180);
				tx+=(speed*updateInterval + acceleration*updateInterval*updateInterval/2)*cos(angle*3.14159/180);
				break;
			case 2://linear motion with sinusoid deviations
				ty=(1+(sin(2*i)/3))*(parentModule()->par("dim_y").doubleValue())/2;
				tx+=(speed*updateInterval + acceleration*updateInterval*updateInterval/2)*cos(angle*3.14159/180);
				break;
			case 3://circle motion
				angle += 4*speed/parentModule()->par("dim_x").doubleValue() * updateInterval;
				tx = parentModule()->par("dim_x").doubleValue()*(2+cos(angle))/4;
				ty = parentModule()->par("dim_x").doubleValue()*(2+sin(angle))/4;
				break;
			}
			speed+=acceleration*updateInterval;
			if (flag){
				appry=C[0]+C[1]*tx+C[2]*tx*tx+C[3]*tx*tx*tx;
				/*trajectoryX.record(tx);
				trajectoryY.record(ty);
				trajectoryapprY.record(appry);*/
				//ev<<"x="<<tx<<"| y_th="<<ty<<" y_ap"<<appry<<"\n";
				if(max_dy<abs(appry-ty)/ty) max_dy=abs(appry-ty)/ty;
				avg_err+=abs(appry-ty)/ty;
				avg_msr++;
			}
		}
		/*if (flag) {
			ev<<"Current MAX_dy="<<max_dy<<" appry="<<appry<<" ty="<<ty<<" tx="<<tx<<"\n";
			ev<<"Coeff: C[0]="<<C[0]<<" C[1]="<<C[1]<<" C[2]="<<C[2]<<" C[3]="<<C[3]<<"\n";
		}*/
	return max_dy;
}
void Node::finish()
{ 
	/*recordScalar("#sent", numSent);
	recordScalar("#target_mes_received", numTReceived);
	recordScalar("#node_mes_received", numNReceived);
	recordScalar("#calculations", numCalc);*/
	char display[30];
	sprintf(display, "%d", numCalc);
	displayString().setTagArg("t",0,display);
	data.clear();
	//Let MASTERNODE to draw errors table
	if (index()==masterNodeID){
		unsigned char counter=0;
		ev<<"Track INFO & ERRORS:\n";
		double max_dx=0,max_dy=0,time0,time1=0,max_err=0;
		ev<<"Time0 "<<m_track.begin()->first<<"\n";//----------<--------DELETE!
		char out[120];
		double x[appr_points],y[appr_points];
		avg_err=0;avg_msr=0;

		double tx0=tx,ty0=ty,speed0=speed,angle0=angle;
		for (trackiter i=m_track.begin();i!=--m_track.end();i++){
			x[counter]=i->second.x;/*Taking appr_points of (x,y) pairs for local approximation*/
			y[counter]=i->second.y;
			max_dx=max(max_dx,i->second.err_x);max_dy=max(max_dy,i->second.err_y);
			sprintf(out,"x=%e | y=%e | err_x=%e | err_y=%e\n",i->second.x,i->second.y,i->second.err_x,i->second.err_y);
			//ev<<out;
			counter++;
			switch (counter){
				case 1: time0=i->first;
					calculateApproximationError(time1,time0,max_err,false);break;
				case appr_points: time1=i->first;//ev<<"Time0: "<<time0<<"\n";ev<<"Time1: "<<time1<<"\n";
					counter=0;
					if(!setApproximationCoefficients(x,y)) max_err=calculateApproximationError(time0,time1,max_err,true);
						else calculateApproximationError(time0,time1,max_err,false);//else is done when we can't calculate approximation
			}
		}
		
		ev<<"MAX dx="<<max_dx<<"% MAX dy="<<max_dy<<"%\n";
		ev<<"-------------\nPolynomial approximation:\n";
		ev<<"MAX dy="<<100*max_err<<"%\n";
		ev<<"AVERAGE dy="<<100*avg_err/avg_msr<<"%\n-------------------\n";

		//counter=0;max_err=0;avg_err=0;avg_msr=0;
		//tx=tx0;ty=ty0;speed=speed0;angle=angle0;
		//time0=0;time1=0;
		//for (trackiter i=m_track.begin();i!=--m_track.end();i++){
		//	x[counter]=i->second.x;/*Taking appr_points of (x,y) pairs for local approximation*/
		//	y[counter]=i->second.y;
		//	counter++;
		//	switch (counter){
		//		case 1: time0=i->first;
		//			calculateApproximationError(time1,time0,0,false);
		//			break;
		//		case 2: time1=i->first;
		//			counter=0;
		//			if(setLinearInterpolationCoefficients(x,y)) max_err=calculateApproximationError(time0,time1,max_err,true);
		//			else {m_track.clear();goto exit; calculateApproximationError(time0,time1,0,false);}//else is done when we can't calculate approximation
		//	}			

		//}
		//ev<<"Simtime "<<simTime()<<"  j="<<m_track.size()-1<<"\n";
		//ev<<"-------------\nLinear interpolation:\n";
		//ev<<"MAX dy="<<100*max_err<<"%\n";
		//ev<<"AVERAGE dy="<<100*avg_err/avg_msr<<"%\n";
		m_track.clear();
	}
	if(false){exit: ev<<"\n----------------------------------------\nATTENTION:\nTime resolution is TOO HIGH!!!\n breaking linear interpolation calculation...\n--------------------------------";}
	/*Now we simply remove all connections so we can see
	results, drawn above each node of sensor network*/
	for (int i=0; i<gate("out")->size(); i++) gate("out",i)->disconnect();
}