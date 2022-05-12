#include "read_graph_data.h"

#include "graph_optimize.h"

#include <vector>
#include <map>

int main()
{
	std::cout<<"------------------ POSE GRAPH TEST ------------------"<<std::endl;
	simulation::ReadGraphData simu;
	simu.openEdgeFile( "./data/killian-e.dat" );
	simu.openVertexFile( "./data/killian-v.dat" );

	std::vector<int> vertex_ids;
	std::vector<Eigen::Vector3f> vertex_poses;

	std::vector<int> edge_from_ids;
	std::vector<int> edge_to_ids;
	std::vector<Eigen::Vector3f> edge_means;

	while( !simu.endOfEdgeFile() ){
//		std::cout<<"edge frame count : "<<simu.getEdgeCount()<<std::endl;
		int id_from = 0;
		int id_to = 0;
		Eigen::Vector3f mean( 0.0, 0.0, 0.0 );
		
		simu.readAEdge( id_from, id_to, mean );
		
		edge_from_ids.push_back( id_from );	
		edge_to_ids.push_back( id_to );
		edge_means.push_back( mean );
	}

	while( !simu.endOfVertexFile() ){
//		std::cout<<"vertex frame count : "<<simu.getVertexCount()<<std::endl;
		int vertex_id = 0;
		Eigen::Vector3f pose( 0.0, 0.0, 0.0 );
		
		simu.readAVertex( vertex_id, pose );
		
		vertex_ids.push_back( vertex_id );
		vertex_poses.push_back( pose );
		
	}
	
	// test
	for( int i = 0; i < vertex_ids.size(); i ++ ) {
		std::cout<<"vertex: "<<vertex_ids[i]<<", pose: "<<std::endl<<vertex_poses[i]<<std::endl;
	}

	for( int i = 0; i < edge_from_ids.size(); i ++ ){
		std::cout<<"id from : "<<edge_from_ids[i]<<" to : "<<edge_to_ids[i]<<std::endl;
	}

	return 0;
}
