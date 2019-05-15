//  Created by Irvanda Kurniadi Virdaus
//  Copyright (c) 2016 Irvanda Kurniadi Virdaus. All rights reserved.
//

// Broadcast simulator for single line road


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define max_dist 5000
#define n_dist 40
#define n_msg 4	// number of message generated 1, 2, 4, 8, 16
#define n_lane 4
#define iter 100

#define CWmin	32
#define CWmax	1024

struct UE { // UE = User Entity / node
    //int cwValue; // cwValue = contention window value
    // int cwSize; // cwSize = contention window size: 8, 16, 32, 64, 128, 256
	int x_pos; // position x
	int y_pos; // position y
};

struct Packet { // packet 
	int nslot; // number of slot
	int nhop; // number of hops
	int col;
	int p_id; // packet id
	int cwValue[(max_dist/n_dist)*n_lane]; // selected CW
	int drop; // packet arrived at last destination
	bool end; // end of packet
	bool is_hear[(max_dist/n_dist)*n_lane];
	bool is_transmit;
	
	bool is_source; // is_source node 
	bool is_sent; // try
	bool is_drop;
	bool transmit[(max_dist/n_dist)*n_lane]; // if the node is transmitting, transmit = true, else false
};

int compute_cw_gd(double dist, double distmax){
	int tmp;
	if (dist < distmax)
		tmp = (int) ((((distmax-dist)/distmax)*(CWmax-CWmin))+CWmin);
	else
		tmp = (int) ((((dist-distmax)/distmax)*(CWmax-CWmin))+CWmin);
	return tmp;
}

int compute_cw_rss(double Pr, double Pr_min, double RxThres){
	int tmp;
	double range = Pr_min - RxThres;
	if (Pr < RxThres)
		tmp = (int) ((((RxThres-Pr)/range)*(CWmax-CWmin))+CWmin);
	else
		tmp = (int) ((((Pr-RxThres)/range)*(CWmax-CWmin))+CWmin);
	return tmp;
}

double compute_distance(int x, int y, int x_s, int y_s, bool is_abs){
	double distance = sqrt(pow((x-x_s),2) + pow((y-y_s),2));
	if (x < x_s && !is_abs)
		distance = -distance;
	return distance;
}

double Friis(double Pt, double Gt, double Gr, double lambda, double L, double d)
{
        /*
         * Friis free space propagation equation:
         *
         *       Pt * Gt * Gr * (lambda^2)
         *   P = --------------------------
         *       (4 *pi * d)^2 * L
         */
  double M = lambda / (4 * M_PI * d);
  return (Pt * Gt * Gr * (M * M)) / L;
}

double TwoRay(double Pt, double Gt, double Gr, double ht, double hr, double L, double d, double lambda)
{
        /*
         *  if d < crossover_dist, use Friis free space model
         *  if d >= crossover_dist, use two ray model
         *
         *  Two-ray ground reflection model.
         *
         *	     Pt * Gt * Gr * (ht^2 * hr^2)
         *  Pr = ----------------------------
         *           d^4 * L
         *
         * The original equation in Rappaport's book assumes L = 1.
         * To be consistant with the free space equation, L is added here.
         */

	double Pr;  // received power
	double crossover_dist = (4 * M_PI * ht * hr) / lambda;

	if (d < crossover_dist)
		Pr = Friis(Pt, Gt, Gr, lambda, L, d);
	else
		Pr = Pt * Gt * Gr * (hr * hr * ht * ht) / (d * d * d * d * L);
		
	return Pr;
}


int main()
{
	//Printout to file
	FILE* nhopFile = fopen ("nhopFile_low_4.tr", "a+");
	FILE* collisionFile = fopen("collisionFile_low_4.tr", "a+");
	FILE* propDelayFile = fopen("nslotFile_low_4.tr","a+");

    // Global variables
	int b_type; // type of broadcast -> 0: simple CW = 32, 1: simple CW = 1024, 2: gd, 3: rss, 4: lab 25, 5: lab 50  
	//int b_type	= 0; // 0: convMin, 1: convMax, 2: Distance, 3: RSS, 4: Area 25, 5: area 50
	int r_node[2]; // relay node position
	double node_dist;
	
	int t_range = 250; // transmission range
	int t_msg =  80; // 10 slots ~ 2 Mbps for 200 bytes = 100 us
	int t_difs =  6 ; // 6 slots ~ 60 us around 58 us
	
    int n = (max_dist/n_dist) * n_lane; // n is number of users in the network
    
    int simSlot = 100000; // Number of time slots in the simulation / duration of simulation (1 s)
                            // 1 slot = 10 us ~ 13 us
                            // 1 seconds = 100,000 slots
	int end_slot;
	
    //GD variables
	double gd_dist;
	
	//LAB variables
	double rad;
	double D_dist;
	double px;
	double py;
	double lab_dist;
	
	//RSS variables
	double Pt = 0.28183815;            // transmit power
	double Gt = 1.0;               // transmit antenna gain
	double Gr = 1.0;               // receive antenna
	double freq = 914.0e6;         // frequency
	double sysLoss = 1.0;  
	double lambda = 3.0e8/freq;
	// for two-ray model
	double ht = 1.5;               // transmit antenna height
	double hr = 1.5;               // receive antenna height
	double Pr; // receiving power
	double RxThres = TwoRay(Pt, Gt, Gr, ht, hr, sysLoss, (double)t_range, lambda);
	double Pr_min = TwoRay(Pt, Gt, Gr, ht, hr, sysLoss,(double)2.0, lambda);
	
    struct UE node[n];
	struct Packet p[n_msg];
	
	int i, j, slot, msg, index, cwSize;
	
	printf("num_nodes: %d \n",n);
	
	//initialize nodes
	for (i = 0; i < n/n_lane; i++){
		for (j = 0;j < n_lane; j++){
			int n_id = j + (i * n_lane);
			node[n_id].x_pos = (i*n_dist)%max_dist+n_dist;
			node[n_id].y_pos = j * 2; //single line
		}
	}
	
	int seed = 2342;
   
	// SIMULATION
	int t_nhop[n_msg] = {};
	int t_nslot[n_msg] = {};
	int nhop_iter[n_msg*iter][6] = {};
	int nslot_iter[n_msg*iter][6] = {};
	int col_iter[n_msg*iter][6] = {};
	int n_receive[n_msg];
	int rem_msg; // remaining msg after dropping
		
	int noTxNodes;
    
	for (b_type = 0; b_type < 6; b_type++){
	
		for (msg = 0; msg < n_msg; msg++){
			
				p[msg].drop = 0;
		}
		
		if (b_type == 4){
			
			rad = 25.0;
			D_dist = 225.0;
			
		}
		else if(b_type == 5){
			
			rad = 50.0;
			D_dist = 200.0;
			
		}

		for (i = 0; i < iter; i++)
		{
			
			srand(seed);
			
			seed = seed + 11;
			
			bool end = false;
			
			int collision = 0;
			
			rem_msg = n_msg;
			
			
			for (j = 0; j < n_msg; j++){
				p[j].nslot = 0;
				p[j].nhop = 0;
				p[j].p_id = j;
				p[j].is_source = true;
				p[j].cwValue[n];
				p[j].is_drop = false;
				p[j].is_transmit = false;
				p[j].end = false;
				p[j].col = 0;
				for (index = 0; index < n; index++){
					
					p[j].is_hear[index] = false;
				}
				
				// printf("nslot : %d, nhop : %d, p_id: %d, is_source: %d \n",p[j].nslot, p[j].nhop, p[j].p_id, p[j].is_source); 
			}
			
			
			for (slot = 0; slot <= simSlot; slot++){
						
				// Set CW value for each nodes
				
				for (msg = 0; msg < n_msg; msg++){
					
					if (p[msg].is_source){
						
						r_node[0] = 0; //set the relay position to 0
						r_node[1] = 0;
						p[msg].is_sent = true;
						p[msg].is_source = false;
						
						slot = slot + t_msg + t_difs;
						p[msg].nslot = p[msg].nslot + t_msg + t_difs;
						// printf("check 1 \n");
					}
					
					if (p[msg].is_sent){
						
						n_receive[msg] = 0;
						
						px = r_node[0] + D_dist;
						py = r_node[1];
						
						p[msg].nhop = p[msg].nhop + 1;
												
						for (index = 0; index < n; index++) {
							
							//set the CW values only in the transmission range to select the relay node
							
							node_dist = compute_distance(node[index].x_pos, node[index].y_pos, r_node[0], r_node[1], false);
							
												
							if (node_dist <= t_range && node_dist > 0 && !p[msg].is_hear[index]){
								
								if (b_type == 0){
									cwSize = (CWmin + 32)-1;
									p[msg].cwValue[index] = rand() % cwSize;
									n_receive[msg]++;
								}
								else if(b_type == 1){
									cwSize = CWmax-1;
									p[msg].cwValue[index] = rand() % cwSize;
									n_receive[msg]++;
								}
								else if(b_type == 2){
									cwSize = compute_cw_gd(node_dist, t_range);
									p[msg].cwValue[index] = rand() % cwSize;
									p[msg].is_hear[index] = true;
									n_receive[msg]++;
								}
								else if(b_type == 3){
									Pr = TwoRay(Pt, Gt, Gr, ht, hr, sysLoss, node_dist, lambda);
									
									cwSize = compute_cw_rss(Pr, Pr_min, RxThres);
									p[msg].cwValue[index] = rand() % cwSize;
									p[msg].is_hear[index] = true;
									n_receive[msg]++;
								}
								else {
									cwSize = (CWmin + 32) - 1;
									lab_dist = compute_distance(px, py, node[index].x_pos, node[index].y_pos, true);
									
									// sqrt(pow((node[index].x_pos-px),2) + pow((node[index].y_pos-py),2));
									// printf("lab_dist = %.2f \n", lab_dist);
									if (lab_dist < rad){
										
										p[msg].cwValue[index] = rand() % cwSize;
										n_receive[msg]++;
										
									}
									else{
										
										p[msg].cwValue[index] = -1; // -1 meand do not rebroadcast. not relay node
									}
									p[msg].is_hear[index] = true;
										
								}
							
								// p[msg].nslot = p[msg].nslot + p[msg].cwValue[index];
								p[msg].transmit[index] = false;
								p[msg].is_transmit = false;
								p[msg].is_sent = false;
								
								//debug
								// printf("r, msg_id %d, node %d, CWvalue %d, time_slot %d \n", msg, index, p[msg].cwValue[index], slot);
								
														
								if (index == n-1){
									p[msg].end = true;
									rem_msg--;
									
								}
								
								// printf("r, msg_id %d, node %d, iter %d, rem_msg %d, CWvalue %d, n_receive %d \n", msg, index, i, rem_msg, p[msg].cwValue[index], n_receive[msg]);
								
								if (rem_msg == 0){
									end = true;
									break;
								}
																		
								// printf("check 2 \n");
							}
						
							else {
								//set CW = -1 for other nodes
								p[msg].cwValue[index] = -1; // -1 means will not be transmitted
								p[msg].transmit[index] = false;
								
							}
							
							
							// printf("check 3 \n");
						}
						
						// printf("iter %d b_type %d msg %d n_receive %d \n", i, b_type, msg, n_receive[msg]);
					}
				} // end of for msg
				
				if (end) //end the simulation if reaching the last node
					break;
					
				//transmission
			
				// Check inside the array NODE if there is any 0 values
			
				noTxNodes = 0;
				
				for (index = 0; index < n; index++) {
				
					for (msg = 0; msg < n_msg; msg++){
					
						if (p[msg].cwValue[index] == 0 && !p[msg].end)
						{
							p[msg].transmit[index] = true;
							noTxNodes = noTxNodes + 1;
							p[msg].is_transmit = true;
							// printf("detect CW = 0 and noTxNodes = %d \n", noTxNodes);
							
							n_receive[msg] = n_receive[msg] - 1;
							
							// printf("check n_receive = %d msg %d node %d noTxNodes %d in check 2 \n", n_receive[msg], msg, index, noTxNodes);
													
						}
											
					}
					
				}
				
				for (msg = 0; msg < n_msg; msg++){
						
					if (n_receive[msg] == 0 && noTxNodes > 1 && !p[msg].is_drop && !p[msg].end){
						// printf("Packet %d DROPED \n", msg);
						p[msg].drop++;
						p[msg].is_drop = true;
						rem_msg--;
						end = true;
						
						// printf("Remaining message in collision checking is %d dropped with n_receive[%d] = %d\n", rem_msg, msg, n_receive[msg]);
					}
					
					if (rem_msg == 0)
						slot = simSlot + 1;

				}
				
				if(end)
					break;
				
				// In case of collision occurs:
				if (noTxNodes > 1){
				
					collision++;
					slot = slot + t_msg + t_difs;
					
					for (msg = 0; msg < n_msg; msg++){
					
							p[msg].nslot = p[msg].nslot + t_msg + t_difs;
							p[msg].is_transmit = false;
							
							if(p[msg].is_transmit)
								p[msg].col = p[msg].col + 1;
					}
				
					// printf("collision occur \n");
					
					for (index = 0; index < n; index++) {
						for (msg = 0; msg < n_msg; msg++){
							if (p[msg].transmit[index])
							{
								p[msg].transmit[index] = false;
								// collision++;
								// slot = slot + t_msg + t_difs;
								// p[msg].nslot = p[msg].nslot + t_msg + t_difs;
								p[msg].cwValue[index] = -1; // -1 means packet dropped
								
							}
						}
					}
				}
				
				else if (noTxNodes == 1)
					
				{
					
					// printf("Successful transmission \n");
					// transmission is successful
					for (index = 0; index < n; index++) {
						for (msg = 0; msg < n_msg; msg++){
							if (p[msg].transmit[index])
							{
								p[msg].transmit[index] = false;
								// success++;
								slot = slot + t_msg + t_difs;
								p[msg].nslot = p[msg].nslot + t_msg + t_difs;
								
								r_node[0] = node[index].x_pos;
								r_node[1] = node[index].y_pos;
								p[msg].is_sent = true;
								p[msg].cwValue[index] = -1; // -1 means packet transmitted
								
								//debug
								// printf("f, msg_id %d, node %d, time_Slot %d, iter %d \n", msg, index, slot, i);
								
							}
						}
					}		
				
				}
				
				else if (noTxNodes == 0)
					
				{
				
					for (msg = 0; msg < n_msg; msg++){
						p[msg].nslot = p[msg].nslot + 1;
					}
					// printf("No transmission - backoff timer \n");
					for (index = 0; index < n; index++) {
						
						for (msg = 0; msg < n_msg; msg++){
					
							if (p[msg].cwValue[index] > 0)
							
								p[msg].cwValue[index] = p[msg].cwValue[index] - 1;
								
						
						}
					}
				}
				
				 //end of 1 slot simulation
			
			} // end of simulation time
			for (msg = 0; msg < n_msg; msg++){
				// if (!p[msg].is_drop){
					nhop_iter[(msg*iter)+i][b_type] = p[msg].nhop;
					nslot_iter[(msg*iter)+i][b_type] = p[msg].nslot;
					col_iter[(msg*iter)+i][b_type] = p[msg].col;
					// col_iter[(msg*iter)+i][b_type] = collision;
				// }
				// else{
					// nhop_iter[(msg*iter)+i][b_type] = 0;
					// nslot_iter[(msg*iter)+i][b_type] = 0;
					// col_iter[(msg*iter)+i][b_type] = 0;
				// }
			}
			
		} // end of iteration
		
		
		double avg_rcv = 0.0;
		for (msg = 0; msg < n_msg; msg++){
			
			avg_rcv = avg_rcv + (double)(100-p[msg].drop)/iter;
		
			// printf("b_type %d msg %d Reception rate =  %.2f \n", b_type, msg, (double)(100-p[msg].drop)/iter);
			
		}
		printf("b_type %d Reception rate %.2f \n", b_type, avg_rcv/n_msg);
			
	} // end of b_type
	for (msg = 0; msg < n_msg; msg++){
	
		for (i = 0; i < iter; i++){
			
			// fprintf (nhopFile,"msg%d iter%d ", msg, i);
			// fprintf (propDelayFile,"msg%d iter%d ", msg, i);
			// fprintf (collisionFile,"msg%d iter%d ", msg, i);
			for (j = 0; j < 6; j++){
			
				fprintf (nhopFile,"%d ",nhop_iter[(msg*iter)+i][j]);
				fprintf (propDelayFile,"%d ", nslot_iter[(msg*iter)+i][j]);
				fprintf (collisionFile,"%d ", col_iter[(msg*iter)+i][j]);
				
			}
			
			fprintf (nhopFile,"\n");
			fprintf (propDelayFile,"\n");
			fprintf (collisionFile,"\n");
			
		}
	}


	fclose(nhopFile);
	fclose(propDelayFile);
	fclose(collisionFile);
    return 0;
} // end of main