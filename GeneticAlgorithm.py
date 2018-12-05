#!/usr/bin/env python
# license removed for brevity
import rospy
import array
from std_msgs.msg import String
from visualization_msgs.msg import *
import numpy
import matplotlib.pyplot as plt
from matplotlib import path
from math import *
import copy

from random import randint
import math

class Farm:
    
    def __init__(self, water):
        self.uid = uid
        self.x = x
        self.y = y

class Drone:

    def __init__(self,uid,fuel,water):
        self.uid = uid 
        self.fuel = fuel
        self.water = water

class Gene:
    
    def __init__(self, drone, nodes_to_visit):
        self.drone = drone
        self.nodes_to_visit = nodes_to_visit

    def __str__(self):
        return "x: %s\ny: %s\nNodeList: %s\n" % (self.x, self.y, self.NodeList)


class Chromosome:

    def __init__(self, genes, FarmCount):
        self.set = genes
        self.FarmCount = FarmCount
        self.FitnessValue = 0

    def Copy(self):
        newSet = []
        for i in range(len(self.set)):
            newSet.append(self.set[i])
        return Chromosome(newSet, self.FarmCount)
    
    def __lt__(self, other):
        return self.FitnessValue < other.FitnessValue
    
    def __eq__(self, other):
        return self.FitnessValue == other.FitnessValue

    # returns a new Chromosome
    def Mutation(self, chromosome_set, temp):
        a = randint(0,499)    
        chromosome1 = chromosome_set[a]
        b = randint(0,3)
        c = randint(0,3)
        while b == c:
            c = randint(0,3)
        
        gene1 = chromosome1[b]
        gene2 = chromosome1[c]
        d = randint(0,len(gene1.nodes_to_visit))
        e = randint(0,len(gene2.nodes_to_visit))
        t = gene1.nodes_to_visit[d]
        gene1.nodes_to_visit[d] = gene2.nodes_to_visit[e]
        gene2.nodes_to_visit[e] = t

        new_chromosome = Chromosome(chromosome1)
        temp.append(new_chromosome)

    def Crossover(self, other):
        newGenes1 = []
        newGenes2 = []
        for i in range(len(self.Genes)):
            newGenes1.append(Gene(self.Genes[i].Course, self.Genes[i].Professor, other.Genes[i].Room, other.Genes[i].Day, other.Genes[i].Slot))
            newGenes2.append(Gene(other.Genes[i].Course, other.Genes[i].Professor, self.Genes[i].Room, self.Genes[i].Day, self.Genes[i].Slot))
        ch1 = Chromosome(newGenes1, self.RoomCount, self.DaysCount, self.TimeSlotsCount)
        ch2 = Chromosome(newGenes2, self.RoomCount, self.DaysCount, self.TimeSlotsCount)
        return [ch1, ch2]

    def FitnessFunction(self):
        Arr = self.Genes.nodes_to_visit 
        index = 0           
        cost = 0                            
        for i in range(1, len(Arr)):
            cost += FarmList[index][i]
            index = i
        self.FitnessValue = cost

    def PrintChromo(self):
        print(self.Genes)
        print(self.FitnessValue)


class GeneticAlgorithm:

    def __init__(self, BestSelectCount, ChromosomesGenerateCount):
        self.BestSelectCount = BestSelectCount
        self.ChromosomesGenerateCount = ChromosomesGenerateCount
    
    def InitialPopulation(self):
        Population = []
        for i in range(self.PopulationCount):
            arr = []
            for j in range(32):
                if ( j < 16 ):
                    for k in range(3):
                        arr.append(Gene(self.Courses[j], self.Professors[int(self.Courses[j].getProfessorID())], randi(0, 1), randi(0, 4), randi(0, 7)))
                else:
                    for k in range(2):
                        arr.append(Gene(self.Courses[j], self.Professors[int(self.Courses[j].getProfessorID())], randi(0, 1), randi(0, 4), randi(0, 7)))
            Population.append(Chromosome(arr, self.Rooms, self.Days, len(self.TimeSlots)))
        return Population

    def simulate(self):
        Pop = self.InitialPopulation()
        Q = PriorityQueue()

        for i in range(len(Pop)):
            Pop[i].FitnessFunction()
            Q.put((Pop[i].FitnessValue, Pop[i]))

        for i in range(self.IterationCount):
            newPop = []
            value = 0
            for j in range(self.BestSelectCount):
                newPop.append(Q.get()[1])
                value += newPop[j].FitnessValue
            print("FitnessValueChromo(Combined): ",value)
            newGen = []
            # print(len(newPop))
            # if ( i < int(self.IterationCount/2) ):
            #     # CrossOver
            #     v = int(len(newPop)/2)
            #     for j in range(v):
            #         for k in range(self.ChromosomesGenerateCount):
            #             newGen = newGen + newPop[j].Crossover(newPop[j+v])
            #     newPop = newPop + newGen
            # else:
                # Mutation
            for j in range(len(newPop)):
                for k in range(self.ChromosomesGenerateCount):
                    newGen.append(newPop[j].Jumble())
            newPop = newPop + newGen
            Q = PriorityQueue()
            for j in range(len(newPop)):
                newPop[j].FitnessFunction()
                Q.put((newPop[j].FitnessValue, newPop[j]))

        newPop = []

    def evaluate_fitness(self, chromosome_set):
        cost = 0
        for i in range(len(chromosome_set)):
            cost += chromosome_set[i].FitnessFunction()

def moveSheperd( point1,point2,intensity,direction ): 
    #print("point1 : " , point1, "  point2 " , point2, " PPPPPPPPPPPPPPPP")
    point1 = numpy.asarray(point1).reshape(1,2)
    point2 = numpy.asarray(point2).reshape(1,2)   
    if numpy.shape(point1)[0]!=0 and numpy.shape(point2)[0]!=0:
        #difference=numpy.asarray(point1).reshape(1,2)-numpy.asarray(point2).reshape(1,2)
        difference=point1-point2
        #print("point1 : " , point1, "  point2 " , point2, " pppppppppp")   
        norm = numpy.zeros(numpy.shape(point1)[0])
        for i in range(numpy.shape(point1)[0]):
            norm[i] = numpy.sqrt(numpy.square(difference[i][0]) + numpy.square(difference[i][1]))
        l=numpy.shape(point2)
        #print("numpy.shape(point1) " , numpy.shape(point1) , "  ",point1[0], " " , norm , )
        #if numpy.shape(point1)[0]>1:
        for i in range(numpy.shape(point2)[0]):
            #print(i , "  " , " iiiiiiiiiiiiiiiiiiiiiiiiiiiii")
            if norm[i]==0:
                unitvec=[0,0]
            else:
                unitvec=[difference[i][0]/norm[i],difference[i][1]/norm[i]]
            #print("unitvec " , unitvec , " " , )    
            point1[i][0]=point1[i][0]-direction*(intensity)*unitvec[0]
            point2[i][0]=point2[i][0]+direction*(intensity)*unitvec[0]
            point1[i][1]=point1[i][1]-direction*(intensity)*unitvec[1]
            point2[i][1]=point2[i][1]+direction*(intensity)*unitvec[1]      
    return point1,point2

def move( point1,point2,intensity,direction ):    
    if numpy.shape(point1)[0]!=0 and numpy.shape(point2)[0]!=0:
        difference=numpy.asarray(point1)-numpy.asarray(point2)
        #difference=point1-point2
        norm=numpy.linalg.norm(difference,axis=0)
        l=numpy.shape(point2)
        if numpy.ndim(point1)>1:
            for i in range(l[1]):
                if norm[i]==0:
                    unitvec=[0,0]
                else:
                    unitvec=[difference[0][i]/norm[i],difference[1][i]/norm[i]]
                point1[0,i]=point1[0,i]-direction*(intensity)*unitvec[0]
                point2[0,i]=point2[0,i]+direction*(intensity)*unitvec[0]
                point1[1,i]=point1[1,i]-direction*(intensity)*unitvec[1]
                point2[1,i]=point2[1,i]+direction*(intensity)*unitvec[1]
        else:
            if norm==0:
                unitvec=[0,0]
            else:
                unitvec=[difference[0]/norm,difference[1]/norm]
            point1[0]=point1[0]-direction*(intensity)*unitvec[0]
            point2[0]=point2[0]+direction*(intensity)*unitvec[0]
            point1[1]=point1[1]-direction*(intensity)*unitvec[1]
            point2[1]=point2[1]+direction*(intensity)*unitvec[1]
    return point1,point2

def DroneMove(drone,rrep,field,threshold,field_with_drone,drone_on_field,speed,fuel_finished,watered_field):
	#CHECK HOW MANY FIELDS EMPTY:
	for i in range(numpy.shape(field)[1]):
		if(i not in field_with_drone and i not in watered_field):
			dist = numpy.linalg.norm(drone.T - field[:,i],axis = 1)
			occupied_field = numpy.where(dist < threshold)
			if(numpy.size(occupied_field) > 0):
				field_with_drone.append(i)
				drone_on_field.append(occupied_field[0][0])

	count = 0;
	for i in field_with_drone:
		if(i == -1):
			count = count + 1

	#Make new set of fields for unused ones:
	if(len(field_with_drone) > 0):
		j = 0;
		field_new = numpy.zeros([2,numpy.shape(field)[1] - len(field_with_drone)])
		for i in range(numpy.shape(field)[1]):
			if(i not in field_with_drone):
				field_new[:,j] = field[:,i];
				j = j+1;
	else:
		field_new = field			

	#print("field_new: " , field_new)-
	print("field_with_drone " , field_with_drone)
	print("drone_on_field: " , drone_on_field)	
	#print("old: " , drone)
	for i in range(numpy.size(drone[0])):
		if(i not in drone_on_field and numpy.shape(field_new)[1] > 0 and i not in fuel_finished):
			print("drone_on_field")
			dist_dest = numpy.linalg.norm(field_new.T - drone[:,i],axis = 1)
			nearest_dest = numpy.argmin(dist_dest,axis = 0)
			#print("field ", field)
			field_nearest = [field_new[0,nearest_dest],field_new[1,nearest_dest]]
			Sh,a = moveSheperd(numpy.asarray(drone[:,i]).reshape(1,2),field_nearest,speed,1)
			#print("field[nearesst]: " , field[:,nearest_dest])
			dist_drone = numpy.linalg.norm(drone.T - Sh,axis = 1)
			inside_rrep = numpy.where(dist_drone < rrep)
			for j in (inside_rrep[0]):
				#print("Sh: " , Sh)
				S = [drone[0,j],drone[1,j]]
				Sh,a = moveSheperd(Sh,S,speed,-1)

			drone[0,i] = Sh[0][0]
			drone[1,i] = Sh[0][1] 	
	#print("NEW " , drone)	
		

	return drone, field_with_drone, drone_on_field;
	

def drones_return(drone,return_pos,speed):
	for i in range(numpy.size(drone[0])):
		Sh,a = moveSheperd(numpy.asarray(drone[:,i]).reshape(1,2),return_pos,speed,1)
		
		dist_drone = numpy.linalg.norm(drone.T - Sh,axis = 1)
		inside_rrep = numpy.where(dist_drone < 3.0)
		for j in (inside_rrep[0]):
			#print("Sh: " , Sh)
			S = [drone[0,j],drone[1,j]]
			Sh,a = moveSheperd(Sh,S,0.5,-1)		
		drone[0,i] = Sh[0][0]
		drone[1,i] = Sh[0][1]


	return drone	


def main():
	rospy.init_node('EVS', anonymous=True)
	pub_marker = rospy.Publisher('evs_project', MarkerArray , queue_size=1)
	pub_marker_capacity = rospy.Publisher('capacity', MarkerArray , queue_size=1)


	rate = rospy.Rate(10) # 10hz

	no_of_fields = 20;
	

	field_size = [20.0,20.0];
	

	field = numpy.zeros((2,no_of_fields));
	# field[0][0] = -80.0
	# field[1][0] =   0.0

	# field[0][1] = -80.0
	# field[1][1] =  80.0
	
	# field[0][2] =   0.0
	# field[1][2] =   0.0

	# field[0][3] =   0.0
	# field[1][3] =  80.0	

	# field[0][4] =  80.0
	# field[1][4] =   0.0

	# field[0][5] =  80.0
	# field[1][5] =  80.0


	field_Xstart = -160;
	field_Ystart =    0;
	rows = 4; cols = 5; k = 0;
	for i in range(rows):
		for j in range(cols):
			field[0][k] = field_Xstart + j*80;
			field[1][k] = field_Ystart + i*80;	
			k = k + 1;

	#drone = numpy.zeros((2,no_of_drones));
	no_of_drones = 8;
	#drone_size = [3.0,3.0];
	drone_size = 3.0
	Xstart =   0.0;
	Ystart = -120.0;
	dronex = Xstart + numpy.random.rand(1,no_of_drones);
	droney = Ystart + numpy.random.rand(1,no_of_drones);

	drone = numpy.concatenate((dronex,droney),axis=0)
	all_drones = []
	for i in range(no_of_drones):
		all_drones.append(i)
	
	#N = 10; nos = 1; nod = 2

	M = MarkerArray()
	shape = Marker.CUBE;
	
	for i1 in range(no_of_drones + no_of_fields): #to include shepherd too
		marker = Marker()	   
		#print(i1)
		marker.header.frame_id = "/multi_sheep_sheperd";
		marker.header.stamp = rospy.get_rostime();
		

		marker.ns = "fields";
		
		marker.id = i1;
		
		marker.type = shape;
		
		marker.action = Marker.ADD;

		marker.color.a = 1.0
		if(i1 < no_of_fields):
			#print("in the if: " , i1)
			marker.pose.position.x = field[0][i1];            
			marker.pose.position.y = field[1][i1];
			marker.pose.position.z = 0.0;
			marker.scale.x = field_size[0];
			marker.scale.y = field_size[1];
			marker.scale.z = 0.1;
			#marker.color.g = 0.15*(i1+1);
			marker.color.g = 1.0


		elif(i1 >= no_of_fields and i1 < no_of_fields+no_of_drones):		    		    		 
			marker.ns = "drones";
			marker.type = Marker.CYLINDER

			marker.pose.position.x = dronex[0][i1 - no_of_fields];            
			marker.pose.position.y = droney[0][i1 - no_of_fields];
			marker.pose.position.z = 5.0;
			marker.scale.x = drone_size
			marker.scale.y = drone_size
			#marker.color.b = 0.15*(no_of_fields - i1 + 1);
			marker.color.b = 0.0
			marker.color.g = 0.0
			marker.color.r = 1.0

			
		marker.pose.orientation.x = 0.0;
		marker.pose.orientation.y = 0.0;
		marker.pose.orientation.z = 0.0;
		marker.pose.orientation.w = 0.0;
		
		marker.lifetime = rospy.Duration();
		M.markers.append(marker)

	C = MarkerArray()
	capacity_size = 20.0;
	for i1 in range(no_of_drones):
		marker = Marker()	   
		#print(i1)
		marker.header.frame_id = "/multi_sheep_sheperd";
		marker.header.stamp = rospy.get_rostime();
		

		marker.ns = "capacity";
		
		marker.id = i1;
		
		marker.type = marker.CYLINDER;
		
		marker.action = Marker.ADD;

		marker.color.a = 1.0
		marker.color.b = 0.0
		marker.color.g = 1.0
		marker.color.r = 1.0		

		marker.pose.position.x = M.markers[no_of_drones+i1].pose.position.x            
		marker.pose.position.y = M.markers[no_of_drones+i1].pose.position.y
		marker.pose.position.z = 0.0;
		marker.scale.x = capacity_size;
		marker.scale.y = capacity_size;
		marker.scale.z = 0.25;

		marker.pose.orientation.x = 0.0;
		marker.pose.orientation.y = 0.0;
		marker.pose.orientation.z = 0.0;
		marker.pose.orientation.w = 0.0;
		
		marker.lifetime = rospy.Duration();
		C.markers.append(marker)		

	GeneticAlgorithm = GeneticAlgorithm(50, 30)		
		 
	rrep = 20.0;	
	threshold = 10;

	watering_speed = 0.02
	speed = 2.0
	
	fuel = capacity_size
	fuel_for_return_journey = 8.0
	fuel_consumption_rate = 0.08
	
	field_with_drone = []
	drone_on_field = []	
	watered_field = []
	fuel_finished = []

	return_pos = [Xstart,Ystart]

	while not rospy.is_shutdown():
		pub_marker.publish(M)
		pub_marker_capacity.publish(C)
		#print("fields: " , field)
		for i in range(no_of_fields):
			M.markers[i].pose.position.x = field[0][i]
			M.markers[i].pose.position.y = field[1][i]
			#print(i, " fields: " , M.markers[i].pose.position.x , " " ,M.markers[i].pose.position.y)

		
		for i in range(no_of_drones):
			M.markers[no_of_fields + i].pose.position.x = drone[0][i]
			M.markers[no_of_fields + i].pose.position.y = drone[1][i]
			#CAPACITY
			C.markers[i].pose.position.x = drone[0][i]
			C.markers[i].pose.position.y = drone[1][i]
			if(i in drone_on_field): #half the fuel consumed
				new_fuel_consumption_rate = fuel_consumption_rate;
			else:
				new_fuel_consumption_rate = fuel_consumption_rate/2.0;

			#Decreasing only till drone size	
			if(C.markers[i].scale.x > drone_size):
				C.markers[i].scale.x = C.markers[i].scale.x - new_fuel_consumption_rate
				C.markers[i].scale.y = C.markers[i].scale.y - new_fuel_consumption_rate 
			
			if(C.markers[i].scale.x < fuel_for_return_journey and i not in fuel_finished):
				fuel_finished.append(i)			
				

				
			# if(len(fuel_finished) == no_of_drones): #returning 
			# 	C.markers[i].scale.x = C.markers[i].scale.x - new_fuel_consumption_rate
			# 	C.markers[i].scale.y = C.markers[i].scale.y - new_fuel_consumption_rate 				

		fuel_not_finished = numpy.setdiff1d(all_drones,fuel_finished)
		#print("fuel_not_finished: " , fuel_not_finished[1] , "")
		fuel_not_finished_list = []
		#Making a list from arrray so that can pass in the array of drone
		for i in range(numpy.size(fuel_not_finished)):
			fuel_not_finished_list.append(fuel_not_finished[i])

		drones_fuel_not_finished , field_with_drone , drone_on_field = DroneMove(drone[:,fuel_not_finished_list],rrep,field,threshold,field_with_drone ,drone_on_field,speed,
			fuel_finished,watered_field)		
		#DroneMove(drone,rrep,field)
		j = 0;
		for i in fuel_not_finished:
			drone[0,i] = drones_fuel_not_finished[0,j];
			drone[1,i] = drones_fuel_not_finished[1,j];
			j = j+1;

		for i in (field_with_drone):
			ind1 = field_with_drone.index(i)
			ind2 = drone_on_field[ind1] #now if this ind2 is present in fuel_finished, then we will not water
			if(i not in watered_field and  ind2 not in fuel_finished): #watering here
				M.markers[i].color.b = M.markers[i].color.b + watering_speed
				M.markers[i].color.g = M.markers[i].color.g - watering_speed/2.0
 				if(M.markers[i].color.b > 0.9):
					watered_field.append(i)
					ind = field_with_drone.index(i)
					drone_on_field[ind] = -1	
					#if(len(drone_on_field) > 0):
					#	drone_on_field.pop(ind)
		# if(len(fuel_finished) > 0):
		# 	for i in fuel_finished:
		# 		ind = drone_on_field.index(i)
		# 		field_with_drone[ind] = -1


		#CALL DRONES BACK : 
		if(len(watered_field) == no_of_fields or len(fuel_finished) == no_of_drones):
			drone = drones_return(drone,return_pos,speed)
		
		if(len(fuel_finished) > 0):
			print("Entering fuel_finished")
			print(fuel_finished)
			drones_fuel_finished = drones_return(drone[:,fuel_finished],return_pos,speed)
			j = 0;
			for i in fuel_finished:
				drone[0,i] = drones_fuel_finished[0,j];
				drone[1,i] = drones_fuel_finished[1,j];
				j = j+1;
			for i in fuel_finished:
				if(i in drone_on_field):
					drone_on_field[drone_on_field.index(i)] = -1;
					#remove the field from field_with_drone
					#field_with_drone[drone_on_field.index(i)] = -1;
					#field_with_drone = field_with_drone.remove(-1)
					
		#drones_re
		#no_of_drones_returned 	

		#print("yoyo")
		rate.sleep();



if __name__ == '__main__':
	FarmList = [[randint(30,200) for i in range(50)] for j in range(50)]
    farm_capacity = [randint(3,60) for i range(50)]
    no_of_drones = 4
    chromosome_set = []

    for i in range(50):
        for j in range(50):
            if i == j:
                FarmList[i][j] = 0
            else:
                FramList[j][i] = -1
    
    for i in range(500):
        chromosome = []
        cities = []
        count = 0
        while count != 50:
            drone_number = randint(0,3)
            city = randint(0,49)
            if not city in cities:
                count += 1
                cities.append(city)
        
            chromosome_set.append(chromosome)

	try:
		main()
	except rospy.ROSInterruptException:
		pass