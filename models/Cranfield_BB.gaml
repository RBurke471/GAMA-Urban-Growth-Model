/**
 *  vector
 *  Author: administrateur
 *  Description: 
 */

model vector

global {
	//Load in shapefiles for model.
    file road_shapefile <- file("../includes/Roads_Poly.shp");
	file road_line_shapefile <- file("../includes/Roads.shp");
	file building_shapefile <- file("../includes/Future_Buildings.shp");
	file building_block_shapefile <- file ("../includes/Parcels_Only.shp");
	file river_shapefile <- file("../includes/Rivers");
	file bounds_shapefile <- file("../includes/Cranfield_CP.shp");
	file river_line_shapefile <- file("../includes/SP_SurfaceWater_Line_Study_Area.shp");
	graph road_network;
	geometry shape <- envelope(bounds_shapefile); //Create bounding box from geometry.
	float step <- 1 #year;  //Timestep for each simulation cycle.
	int cycle;
	float time;
	bool batch_mode <- false;
	bool save_dist <- false;
	
	float w1 <- 0.2; //Weight value for building density criteria
	float w2 <- 0.0; //
	float w3 <- 1.0; //Weight value for roads and rivers criteria
	
	int n <- 387;
	int nb_build <- 387; //Number of buildings to construct each cycle
	
	float distance_roads_rivers <- 1.0;
	
	
	float crit_roads_rivers_max;
	string dist_path <- "../includes/vector_distances.csv";
	matrix dist;
	
	int year <- 2021; //Start year of simulation.
	date starting_date <- date([year]);
		
	init {
		create road from: road_shapefile;
		create road_line from: road_line_shapefile;
		road_network <- as_edge_graph(road_line);
		create river from: river_shapefile;
		create river_line from: river_line_shapefile;
		create building from: building_shapefile with: [type::string(get("BUILDGTHEM"))];
		create building_block from: building_block_shapefile;
		create Cranfield from: bounds_shapefile;
		if (not save_dist) {
			dist <- matrix(csv_file(dist_path, ";"));
		}
		do building_block_creation;
		if (save_dist) {
			do compute_distances;
		}
	}
	
	reflex global_dynamic when: w1 != 0 or w3 != 0 {
		ask building_block {
			color <- #red;
		}
		list<building_block> bb_to_builds <- building_block where (each.possible_construction) ;
		ask bb_to_builds{
			do compute_criteria;
		}
		crit_roads_rivers_max <- building_block max_of (each.crit_roads_rivers);
			
		ask bb_to_builds{
			do compute_land_value;
		}
		list<building_block> sorted_bb <- shuffle(bb_to_builds) sort_by (- each.constructability); //Creat list of building_block sorted by constructability 
		loop i from: 0 to: min([n, length(sorted_bb)]){
			building_block bb_to_build <- sorted_bb[i];
			ask bb_to_build {
				loop times: nb_build {
					do building_construction;
				}
			}
		}
	}
	
	
	reflex end_simulation when: date = 2100 and not batch_mode{
		do pause;
	}
	 
	action compute_distances { 
		int nb <- length(building_block);
		dist <- 0 as_matrix {nb,nb};
		loop i from: 0 to: nb - 1 {
			building_block bb1 <- building_block[i];
			loop j from: 0 to:  nb - 1   {
				building_block bb2 <- building_block[j];
				if (bb1 = bb2) {
					dist[i,j] <- 0.0;
				} else {
					dist[i,j] <- topology(road_network) distance_between([bb1,bb2]);
				}
		 	}		 	
		 } 
		 save string(dist) to: dist_path type:"text";
		  
	}
	
	action building_block_creation {
		if (not save_dist) {
			ask building_block {
				buildings <- building; //Create layer called buildings from building that overlap the building_block
				do build_empty_space;			}
		}
	}
}

species river {
	aspect geom {
		draw shape color: #blue;
	}
}

species river_line {
	aspect geom {
		draw shape color: #blue;
	}
}
species road {
	aspect geom {
		draw shape color: #yellow;
	}
}

species road_line {
	aspect geom {
		draw shape color: #red;
	}
}

species building {
	int build_year <- 9999; //Give current buildings a value not associated with simulation.
	string type;
	aspect geom {
		draw shape color: #green;
	}
}

species Cranfield {
	string type;
	aspect geom {
		draw shape color: #black;
	}
}




species building_block {
	rgb color <- #red;
	list<building> buildings;
	geometry empty_space; 
	bool possible_construction <- true ;  //boolean to determine if construction is possible
	float free_space_rate;
	float constructability;
	float crit_roads_rivers;
	
	
	action distance_to_roads {
		
	}
	
	action build_empty_space {
		empty_space <- copy(shape); 
		loop bd over: buildings {
			empty_space <- empty_space - (bd + 2.0); //
		}
		
		if (empty_space != nil) {
			list<geometry> geoms_to_keep <- empty_space.geometries where (each != nil and each.area > 200);
			if (not empty(geoms_to_keep)) {
				empty_space <- geometry(geoms_to_keep);
			} else {
				empty_space <- nil;
				possible_construction <- false;
			}
		}
	}
	
	action update_empty_space(building bd) {
		empty_space <- empty_space - (bd + 2.0); //Remove buildings (with 2m buffer) from empty space
		//If empty spaces are not equal to nil then list their geometries greater than 200m2
		if (empty_space != nil) {
			list<geometry> geoms_to_keep <- empty_space.geometries where (each != nil and each.area > 200);
			if (not empty(geoms_to_keep)) {
				empty_space <- geometry(geoms_to_keep);
			} else {
				empty_space <- nil;
				possible_construction <- false;
			}
		}
	}
	
	action compute_criteria {
	
	list roads <- road at_distance distance_roads_rivers;  
	geometry shape_buffer <- shape + distance_roads_rivers;
	crit_roads_rivers <-sum(roads collect (each.shape inter shape_buffer).area) / envelope(shape).area;
		
	}
	//Compute building density, quantity of roads and rivers, and overall constructability.
	action compute_land_value {
		crit_roads_rivers <- crit_roads_rivers / crit_roads_rivers_max;
		constructability <- (crit_roads_rivers * w3)/ (w3); 
	}

	
	
	action building_construction {
		float limit <- world.shape.area; //Limit building construction
		list<building> possible_buildings <- (buildings);
		if empty(possible_buildings) {
			possible_construction<- false;
			write name;
			return;
		}
		bool bd_built <- false;
		building bd <- nil;
		
		loop while: true {
			building one_building <- one_of (possible_buildings where (envelope(each.shape).area < limit));
			if (one_building = nil) {
				break;
			}
			geometry new_building <- copy(one_building.shape);
			float size <- min([new_building.width,new_building.height]) ;
			
			geometry space <- empty_space reduced_by size;
			if ((space != nil) and (space.area > 0.0)) {
				agent closest_road_river <- (road) closest_to space;
				create building with:[ shape:: new_building, location::((closest_road_river closest_points_with space)[1]),type::one_building.type] {
					myself.buildings << self;
					build_year <- current_date.year; //Add simulation year to attribute table (build_year column).
					bd <- self;
					switch type {
					}
				}
				bd_built <- true;
				break;
			} else {
				limit <- envelope(one_building.shape).area;
			}
		}
		
		if (bd_built) {
			do update_empty_space(bd);	
		} else {
			possible_construction<- false;
		}
	}
	
	aspect geom {
		draw shape color:color border: #black;
	}
}



experiment vector type: gui until: ( date = 2050 ) {
	parameter "weight for the quantity of roads and rivers criterion" var:w3 min: 0.0 max: 1.0; 
	
	output {
		display "map" type: opengl ambient_light: 100{
			species Cranfield aspect: geom;
			species building_block aspect: geom;
			species river aspect: geom;
			species river_line aspect: geom;
			species road aspect: geom;
			species building aspect: geom;	
		}
	}
}


experiment Optimization type: batch keep_seed: true repeat: 5 until: ( date = 2010 ) {
	parameter "batch mode" var: batch_mode <- true;
	list<float> vals;
	parameter "weight for the density criterion" var:w1 min: 0.0 max: 1.0 step: 0.1; 
	parameter "weight for the distance to services (education, religion, economic)" var:w2 min: 0.0 max: 1.0 step: 0.1;
	parameter "weight for the quantity of roads and rivers criterion" var:w3 min: 0.0 max: 1.0 step: 0.1; 
	
	method genetic pop_dim: 5 crossover_prob: 0.7 mutation_prob: 0.1 nb_prelim_gen: 1 max_gen: 500  minimize: error ;
	
	reflex save_result {
		vals<< world.error;
		save (string(world.w1) + "," + world.w2 + "," + world.w2 + "," + world.error) to: "E:/GAMA/Workspace/results_vector.csv" type:"text"; 
		if (length(vals) = 5) {
			write "error: " + mean(vals) + " for parameters: w1 = " + w1 + "  w2 = " + w2 + "  w3 = " + w3 ;
			vals <- [];
		}
	}
}

experiment save_distances type: gui {
	parameter "save distances" var: save_dist <- true;
	output {}
}