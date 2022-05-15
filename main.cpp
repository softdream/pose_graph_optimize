#include "read_graph_data.h"

#include "graph_optimize.h"

#include <vector>
#include <map>

#include "utils.h"
#include <chrono>

int main()
{
	std::cout<<"------------------ POSE GRAPH TEST ------------------"<<std::endl;
	simulation::ReadGraphData simu;
	simu.openEdgeFile( "./data/killian-e.dat" );
	simu.openVertexFile( "./data/killian-v.dat" );

	std::vector<int> vertex_ids;
	std::vector<Eigen::Vector3d> vertex_poses;

	std::vector<int> edge_from_ids;
	std::vector<int> edge_to_ids;
	std::vector<Eigen::Vector3d> edge_means;
	std::vector<Eigen::Matrix3d> info_matrixes;

	while( simu.getEdgeCount() < 3995 ){
		//std::cout<<"edge frame count : "<<simu.getEdgeCount()<<std::endl;
		int id_from = 0;
		int id_to = 0;
		Eigen::Vector3d mean( 0.0, 0.0, 0.0 );
		Eigen::Matrix3d info_matrix = Eigen::Matrix3d::Zero();
		

		simu.readAEdge( id_from, id_to, mean, info_matrix );
		
		edge_from_ids.push_back( id_from );	
		edge_to_ids.push_back( id_to );
		edge_means.push_back( mean );
		info_matrixes.push_back( info_matrix );
	}

	while( simu.getVertexCount() < 1941 ){
		//std::cout<<"vertex frame count : "<<simu.getVertexCount()<<std::endl;
		int vertex_id = 0;
		Eigen::Vector3d pose( 0.0, 0.0, 0.0 );
		
		simu.readAVertex( vertex_id, pose );
		
		vertex_ids.push_back( vertex_id );
		vertex_poses.push_back( pose );
		
	}
	
	// test
	/*for( int i = 0; i < vertex_ids.size(); i ++ ) {
		std::cout<<"vertex: "<<vertex_ids[i]<<", pose: "<<std::endl<<vertex_poses[i]<<std::endl;
	}

	for( int i = 0; i < edge_from_ids.size(); i ++ ){
		std::cout<<"id from : "<<edge_from_ids[i]<<" to : "<<edge_to_ids[i]<<std::endl;
	}
	for( int i = 0; i < edge_means.size(); i ++ ){
                std::cout<<"edge_means : "<<i<<std::endl<<edge_means[i]<<std::endl;
        }

	for( int i = 0; i < info_matrixes.size(); i ++ ){
		std::cout<<"info_matrix["<<i<<"] : "<<std::endl<<info_matrixes[i]<<std::endl;
	}
*/
	
	std::cout<<"edge_number : "<<edge_means.size()<<std::endl;
	std::cout<<"vertex number : "<<vertex_ids.size()<<std::endl;

	Utils::displayVertexPoses( vertex_poses, "init" );

	// ------------------- Graph Optimize -------------------- //
	graph::GraphOptimize<double> optimizer( vertex_poses.size() );

	auto beforeTime = std::chrono::steady_clock::now();
	optimizer.execuOptimize( vertex_poses, 
				 edge_from_ids, 
				 edge_to_ids, 
				 edge_means,
				 info_matrixes, 2 );
	auto afterTime = std::chrono::steady_clock::now();
        double duration_millsecond = std::chrono::duration<double, std::milli>(afterTime - beforeTime).count();
        std::cout<<"duration : " << duration_millsecond << "ms" << std::endl;

	std::vector<Eigen::Vector3d> ret_vec = optimizer.getReultVertexPosesVector();

	Utils::displayVertexPoses( ret_vec );
	
	return 0;
}
