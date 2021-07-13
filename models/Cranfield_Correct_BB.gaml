/**
 *  vector
 *  Author: administrateur
 *  Description: 
 */

model vector

global {
    file road_shapefile <- file("../includes/Roads_Poly.shp");
	file road_line_shapefile <- file("../includes/Roads.shp");
	file building_shapefile <- file("../includes/Buildings.shp");
	file building_block_shapefile <- file ("../includes/Building_Parcels_Centre.shp");
	file river_shapefile <- file("../includes/Rivers");
	file bounds_shapefile <- file("../includes/Cranfield_CP.shp");
	graph road_network;
	geometry shape <- envelope(bounds_shapefile); //Create bounding box from geometry.
	float step <- 1 #year;
	list<geometry> old_buildings;
	list<geometry> new_buildings;
	float error;
	list<geometry> new_buildings_init;
	list<geometry> old_buildings_init;
	
	bool batch_mode <- false;
	bool save_dist <- false;
	
	
	float w1 <- 0.2;
	float w2 <- 0.0;
	float w3 <- 1.0;
	
	int n <- 10;
	int nb_build <- 10; //Number of buildings to construct each simulation
	float distance_neighbours <- 50.0;
	float distance_roads_rivers <- 1.0;
	
	
	float crit_density_max;
	float crit_roads_rivers_max;
	string dist_path <- "../includes/vector_distances.csv";
	matrix dist;
	
	init {
		create road from: road_shapefile;
		create road_line from: road_line_shapefile;
		road_network <- as_edge_graph(road_line);
		create river from: river_shapefile;
		create building from: building_shapefile with: [type::string(get("BUILDGTHEM"))];
		create building_block from: building_block_shapefile;
		if (not save_dist) {
			dist <- matrix(csv_file(dist_path, ";"));
		}
		do building_block_creation;
		
	}
	
	reflex global_dynamic when: w1 != 0 or w3 != 0 {
		ask building_block {
			color <- #red;
		}
		list<building_block> bb_to_builds <- building_block where (each.possible_construction); //List building blocks (called bb_to_builds) that are possible to construct in at each cycle.
		ask bb_to_builds{
			do compute_criteria;
		}
		crit_density_max <- building_block max_of (each.crit_density);
		crit_roads_rivers_max <- building_block max_of (each.crit_roads_rivers);
			
		ask bb_to_builds{
			do compute_land_value;
		}
		
		list<building_block> sorted_bb <- shuffle(bb_to_builds) sort_by (each.constructability);  //List sorted building blocks in order of constructability 
		loop i from: 0 to: min([n, length(sorted_bb)]){  //
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
	
	action building_block_creation {
		ask building_block {
			buildings <- building overlapping self;
			do build_empty_space;			
		}
	}
}

species river {
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
	string type;
	aspect geom {
		draw shape color: #green;
	}
}


species building_block {
	rgb color <- #red;
	list<building> buildings;
	geometry empty_space; 
	list<building_block> neigbours_blocks;
	bool possible_construction <- true ; //Boolean to indicate whether construction is possible
	float free_space_rate;
	bool has_education <- false;
	bool has_religious <- false;
	bool has_economic <- false;
	map<building_block, float> distances;
	float constructability;
	float crit_density;
	float crit_roads_rivers;
	
	
	action build_empty_space {
		empty_space <- copy(shape); 
		loop bd over: buildings {
			empty_space <- empty_space - (bd +2.0);
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
		empty_space <- empty_space - (bd + 2.0);
		
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
		crit_density <- empty_space = nil ? 0.0 : (shape.area/empty_space.area);
		
		list roads <- road at_distance distance_roads_rivers;
		list rivers <-river at_distance distance_roads_rivers;
		geometry shape_bufer <- shape + distance_roads_rivers;
		crit_roads_rivers <-(sum(roads collect (each.shape inter shape_bufer).area) + sum(rivers collect (each.shape inter shape_bufer).area)) / envelope(shape).area;
	}
	action compute_land_value {
		crit_density <- crit_density/crit_density_max;
		crit_roads_rivers <- crit_roads_rivers / crit_roads_rivers_max;
		constructability <- (crit_density * w1 + crit_roads_rivers * w3)/ (w1 + w3); //Calculate constructability of the building blocks
	}
	
	
	action building_construction {
		float limit <- world.shape.area;
		list<building> possible_buildings <- (buildings + (neigbours_blocks accumulate each.buildings)) where (each.type = "Residential");
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
					bd <- self;
					switch type {
						match "Education" {myself.has_education <- true;}
						match "Religion" {myself.has_religious <- true;}
						match "Economic activity" {myself.has_economic <- true;}
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



experiment vector type: gui {
	parameter "weight for the density criterion" var:w1 min: 0.0 max: 1.0; 
	parameter "weight for the quantity of roads and rivers criterion" var:w3 min: 0.0 max: 1.0; 
	
	output {
		monitor "Error" value: error;
		display "map" type: opengl ambient_light: 100{
			species building_block aspect: geom;
			species river aspect: geom;
			species road aspect: geom;
			species building aspect: geom;	
		}
	}
}

experiment save_distances type: gui {
	parameter "save distances" var: save_dist <- true;
	output {}
}

