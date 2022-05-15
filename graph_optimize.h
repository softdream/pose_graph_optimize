#ifndef __GRAPH_OPTIMIZE_H
#define __GRAPH_OPTIMIZE_H

#include <vector>
#include <Eigen/Dense>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/Eigen>

#include<Eigen/SparseCholesky>
#include<Eigen/SparseLU>

#include "utils.h"

namespace graph
{

template<typename T>
class GraphOptimize
{
public:
	using DataType = T;
	using Vector3 = typename Eigen::Matrix<DataType, 3, 1>;
	using Matrix3x3 = typename Eigen::Matrix<DataType, 3, 3>;
	using Vector2 = typename Eigen::Matrix<DataType, 2, 1>;
		

	GraphOptimize() = delete;

	explicit GraphOptimize( const int vertex_num )
	{
		H.resize( vertex_num * 3, vertex_num * 3 );
		b.resize( vertex_num * 3, 1 );
		delta_x.resize( vertex_num * 3, 1 );
		x.resize( vertex_num * 3, 1 );
	} 

	~GraphOptimize()
	{

	}

	void setVertexNumber( const int vertex_num )
	{
		H.resize( vertex_num * 3, vertex_num * 3 );
		b.resize( vertex_num * 3, 1 );
		delta_x.resize( vertex_num * 3, 1 );
		x.resize( vertex_num * 3, 1 );
	}

	void execuOptimize( std::vector<Vector3> &v_poses,
                            const std::vector<int> &from_ids,
                            const std::vector<int> &to_ids,
                            const std::vector<Vector3> &e_means,
                            const std::vector<Matrix3x3> &info_matrix,
			    const int max_iterations = 10 )
	{
		for( size_t i = 0; i < v_poses.size(); i ++ ){
			x( i * 3, 0 ) = v_poses[i](0);
			x( i * 3 + 1, 0 ) = v_poses[i](1);
			x( i * 3 + 2, 0 ) = v_poses[i](2);
		}	

		for( int iteration = 0; iteration < max_iterations; iteration ++ ){
			estimateOnce( v_poses, from_ids, to_ids, e_means, info_matrix );
			std::cout<<"estimate iteration : "<<iteration<<std::endl;		
		}

		for( size_t i = 0; i < v_poses.size(); i ++ ){
			//angleNormalize( x( i * 3 + 2 ) );
			DataType theta = x( i * 3 + 2 );
			DataType s = ::sin( theta );
			DataType c = ::cos( theta );
			x( i * 3 + 2 ) = ::atan2( s, c );
		}
	}

	const Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> &getResultVertexPosesMatrix() const
	{
		return x;
	}
	
	const std::vector<Vector3> getReultVertexPosesVector() const
	{
		std::vector<Vector3> ret;
		for( size_t i = 0; i < x.size(); i += 3 ){
			Vector3 pose( x( i ), x( i + 1 ), x( i + 2 ) );
			ret.push_back( pose );
		}

		return ret;
	}
	
private:
	void estimateOnce( std::vector<Vector3> &v_poses,
                           const std::vector<int> &from_ids,
                           const std::vector<int> &to_ids,
                           const std::vector<Vector3> &e_means,
                           const std::vector<Matrix3x3> &info_matrix )
	{
		getHessianDerived( v_poses, from_ids, to_ids, e_means, info_matrix );
		
		// solve the linear system using sparse Cholesky factorization
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<DataType>> solver;
                //Eigen::SparseLU<Eigen::SparseMatrix<DataType>> solver;
                solver.compute( H.sparseView() );

                if (solver.info() != Eigen::Success) {
                        std::cerr << "Decomposition Failed !" << std::endl;
                        return;
                }

                delta_x.setZero();
                delta_x = solver.solve( b );
                if (solver.info() != Eigen::Success) {
                        std::cerr << "Solving Failed !" << std::endl;
                        return;
                }
		
		//std::cout<<"delta_x : "<<std::endl<<delta_x<<std::endl;
                // update
                x += delta_x;
		//std::cout<<"x_new  = "<<std::endl<<x<<std::endl;
		//for( size_t i = 0; i < x.size(); i ++ ){
                //        std::cout<<"x_new[ "<<i<<" ] = "<<x(i, 0)<<std::endl;
                //}

		v_poses.clear();
		for( size_t i = 0; i < x.size(); i += 3 ){
                        Vector3 pose( x( i ), x( i + 1 ), x( i + 2 ) );
                        v_poses.push_back( pose );
                }
	}

	void getHessianDerived( const std::vector<Vector3> &v_poses,
				const std::vector<int> &from_ids,
				const std::vector<int> &to_ids,
				const std::vector<Vector3> &e_means,
				const std::vector<Matrix3x3> &info_matrix )
	{
		H.setZero(); // 1. H <- 0
		b.setZero(); //    b <- 0

		// for all <e_ij, Omega_ij> do:
		for( size_t i = 0; i < e_means.size(); i ++ ){
			//std::cout<<"------------------- edge : "<<i<<" ---------------------"<<std::endl;
			int id_i = from_ids[i];
			int id_j = to_ids[i];
			//std::cout<<"from id_i = "<<id_i<<", to id_j = "<<id_j<<std::endl;
			// compute the Jacobians A_ij and B_ij of the error function
			linearFactors( v_poses, from_ids, to_ids, e_means, i );
			//std::cout<<"caculate linear factors "<<std::endl;	
			// compute the coefficient vector
			b_i = Vector3::Zero();
	                b_j = Vector3::Zero();
        	        H_ii = Matrix3x3::Zero();
                	H_ij = Matrix3x3::Zero();
	                H_ji = Matrix3x3::Zero();
        	        H_jj = Matrix3x3::Zero();
			
			Matrix3x3 omega = info_matrix[i];

			b_i = -A.transpose() * omega * e;
			b_j = -B.transpose() * omega * e;
			// compute the contribution of this constraint to the linear system
			H_ii = A.transpose() * omega * A;
			H_ij = A.transpose() * omega * B;
			//H_ji = B.transpose() * omega * A;
			H_jj = B.transpose() * omega * B;

			if( i == 1499 ){
				std::cout<<"b_i : "<< std::endl<<b_i<<std::endl;
                        	std::cout<<"b_j : "<< std::endl<<b_j<<std::endl;
			}
		
			/*std::cout<<"b_i : "<< std::endl<<b_i<<std::endl;
			std::cout<<"b_j : "<< std::endl<<b_j<<std::endl;
			std::cout<<"H_ii : "<<std::endl<<H_ii<<std::endl;
			std::cout<<"H_ij : "<<std::endl<<H_ij<<std::endl;
			std::cout<<"H_ji : "<<std::endl<<H_ji<<std::endl;
			std::cout<<"H_jj : "<<std::endl<<H_jj<<std::endl;	
			*/
			H.block( id_i * 3, id_i * 3, 3, 3 ) += H_ii;
			H.block( id_i * 3, id_j * 3, 3, 3 ) += H_ij;
			H.block( id_j * 3, id_i * 3, 3, 3 ) += H_ij.transpose();
			H.block( id_j * 3, id_j * 3, 3, 3 ) += H_jj;

			b( id_i * 3, 0 ) += b_i( 0 );
			b( id_i * 3 + 1, 0 ) += b_i( 1 );
			b( id_i * 3 + 2, 0) += b_i( 2 );
		
			b( id_j * 3, 0 ) += b_j( 0 );
			b( id_j * 3 + 1, 0 ) += b_j( 1 );
			b( id_j * 3 + 2, 0 ) += b_j( 2 );

			//std::cout<<"b : "<<b<<std::endl;
			/*for( size_t i = 0; i < b.size(); i ++ ){
                	        std::cout<<"b[ "<<i<<" ] = "<<b(i, 0)<<std::endl;
                	}*/
		}
			
		H( 0, 0 ) = 1;
		H( 1, 1 ) = 1;
		H( 2, 2 ) = 1;
	//	std::cout<<"H all : "<<H<<std::endl;
	//	std::cout<<"b all: "<<b<<std::endl;
	//	for( size_t i = 0; i < b.size(); i ++ ){
	//		std::cout<<"b[ "<<i<<" ] = "<<b(i, 0)<<std::endl;
	//	}
	}


	/** compute the taylor expansion of the error function of the index_th edge
	*   v_poses: vertices position
	*   from_ids, to_ids: edge number
	*   e_means: edges means
	*/
	void linearFactors( const std::vector<Vector3> &v_poses, 
			    const std::vector<int> &from_ids,
			    const std::vector<int> &to_ids,
			    const std::vector<Vector3> &e_means,
			    const int index )
	{
		
		A = Matrix3x3::Zero();
		B = Matrix3x3::Zero();
		e = Vector3::Zero();
	
		int id_i = from_ids[index];
		int id_j = to_ids[index];
		
		Vector3 v_i = v_poses[id_i];
		Vector3 v_j = v_poses[id_j];
		Vector3 z_ij = e_means[index];
		//if( index == 1499 ){
		//	std::cout<<"v_i = "<<id_i<<std::endl<<v_i<<std::endl;
		//	std::cout<<"v_j = "<<id_j<<std::endl<<v_j<<std::endl;	
		//	std::cout<<"z_ij = "<<std::endl<<z_ij<<std::endl;
		//}

		Matrix3x3 zt_ij = v2t( z_ij );
		Matrix3x3 vt_i = v2t( v_i );
		Matrix3x3 vt_j = v2t( v_j );
		//if( index == 1499 ){
		//	std::cout<<"zt_ij = "<<std::endl<<zt_ij<<std::endl;
		//	std::cout<<"vt_i = "<<std::endl<<vt_i<<std::endl;
		//	std::cout<<"vt_j = "<<std::endl<<vt_j<<std::endl;
		//}
		
		Matrix3x3 f_ij = vt_i.inverse() * vt_j;
		//if( index == 1499 ){
		//	std::cout<<"f_ij = "<<std::endl<<f_ij<<std::endl;		
		//}

		DataType theta_i = v_i[2];
		Vector2 t_i( v_i[0], v_i[1] );
		Vector2 t_j( v_j[0], v_j[1] );
	
		Vector2 dt_ij = t_j - t_i;
		//std::cout<<"theta_i = "<<theta_i<<", dt_ij = "<<std::endl<<dt_ij<<std::endl;

		DataType si = ::sin( theta_i );
		DataType ci = ::cos( theta_i );

		A << -ci, -si, -si * dt_ij[0] + ci * dt_ij[1],
		      si, -ci, -ci * dt_ij[0] - si * dt_ij[1],
		       0,   0,  -1;
	
		B << ci, si, 0,
		    -si, ci, 0,
		      0,  0, 1;
		//if( index == 1499 ){
		//	std::cout<<"A_pre = "<<std::endl<<A<<std::endl;
		//	std::cout<<"B_pre = "<<std::endl<<B<<std::endl;
		//}	
		Matrix3x3 zt_ij_inv = zt_ij.inverse();
		//std::cout<<"zt_ij_inv = "<<std::endl<<zt_ij_inv<<std::endl;
	
		e = t2v( zt_ij_inv * f_ij );
		//if( index == 1499 ){
		//	std::cout<<"e = "<<e.transpose()<<std::endl;	
		//}	
	
		zt_ij_inv( 0, 2 ) = 0;
		zt_ij_inv( 1, 2 ) = 0;
	
		A = zt_ij_inv * A;
		B = zt_ij_inv * B;
		//if( index == 1499 ){
		//	std::cout<<"A = "<<std::endl<<A<<std::endl;
		//	std::cout<<"B = "<<std::endl<<B<<std::endl;
		//}
	}

private:
	const Matrix3x3 v2t( const Vector3 &v )
	{
		Matrix3x3 A;
	
		A << ::cos( v[2] ), -::sin( v[2] ), v[0],
		     ::sin( v[2] ),  ::cos( v[2] ), v[1],
			0,                 0,        1;

		return A;
	}

	const Vector3 t2v( const Matrix3x3 &A )
	{
		Vector3 v;
		
		v[0] = A( 0, 2 );
		v[1] = A( 1, 2 );
		v[2] = ::atan2( A( 1, 0 ), A( 0, 0 ) );
	
		return v;
	}

	void angleNormalize( DataType &angle )
	{
        	if( angle >= M_PI ) {
               		angle -= 2 * M_PI;
        	}

        	if( angle <= -M_PI ){
                	angle += 2 * M_PI;
        	}
	}

private:
	//Eigen::SparseMatrix<DataType> H;	
	Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> H;	

	Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> b;
	Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> delta_x;
	Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> x;

	Matrix3x3 A;
	Matrix3x3 B;
	Vector3 e;

	Vector3 b_i;
	Vector3 b_j;
	Matrix3x3 H_ii;
	Matrix3x3 H_ij;
	Matrix3x3 H_ji;
	Matrix3x3 H_jj; 
};

}

#endif
