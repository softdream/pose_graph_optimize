#ifndef __GRAPH_OPTIMIZE_H
#define __GRAPH_OPTIMIZE_H

#include <vector>
#include <Eigen/Dense>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/Eigen>

#include<Eigen/SparseCholesky>

#define PI 3.141592653

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

	void execuOptimize( const std::vector<Vector3> &v_poses,
                            const std::vector<int> &from_ids,
                            const std::vector<int> &to_ids,
                            const std::vector<Vector3> &e_means,
                            const Matrix3x3 &info_matrix,
			    const int max_iterations = 10 )
	{
		for( size_t i = 0; i < v_poses.size(); i ++ ){
			x( i * 3, 0 ) = v_poses[i](0);
			x( i * 3 + 1, 0 ) = v_poses[i](1);
			x( i * 3 + 2, 0 ) = v_poses[i](2);
		}	

		for( int iteration = 0; iteration < max_iterations; iteration ++ ){
			estimateOnce( v_poses, from_ids, to_ids, e_means, info_matrix );
		}

		for( size_t i = 0; i < v_poses.size(); i ++ ){
			angleNormalize( x( i * 3 + 2, 0 ) );
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
	void estimateOnce( const std::vector<Vector3> &v_poses,
                           const std::vector<int> &from_ids,
                           const std::vector<int> &to_ids,
                           const std::vector<Vector3> &e_means,
                           const Matrix3x3 &info_matrix )
	{
		getHessianDerived( v_poses, from_ids, to_ids, e_means, info_matrix );
		
		// solve the linear system using sparse Cholesky factorization
		Eigen::SimplicialLLT<Eigen::SparseMatrix<DataType>> solver;
		solver.compute( H );
		
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
		
		// update 
		x += delta_x;
	}

	void getHessianDerived( const std::vector<Vector3> &v_poses,
				const std::vector<int> &from_ids,
				const std::vector<int> &to_ids,
				const std::vector<Vector3> &e_means,
				const Matrix3x3 &info_matrix )
	{
		//Eigen::SparseMatrix<DataType> H( v_poses.size() * 3, v_poses.size() * 3 );
		b_i = Vector3::Zero();
		b_j = Vector3::Zero();
		H_ii = Matrix3x3::Zero();
		H_ij = Matrix3x3::Zero();
		H_ji = Matrix3x3::Zero();
		H_jj = Matrix3x3::Zero();		

		std::vector <Eigen::Triplet<DataType>> H_triplets;		

		H.setZero(); // 1. H <- 0
		b.setZero(); //    b <- 0

		// for all <e_ij, Omega_ij> do:
		for( size_t i = 0; i < e_means.size(); i ++ ){
			int id_i = from_ids[i];
			int id_j = to_ids[i];

			// compute the Jacobians A_ij and B_ij of the error function
			linearFactors( v_poses, from_ids, to_ids, e_means, i );
			
			// compute the coefficient vector
			b_i = -A.transpose() * info_matrix * e;
			b_j = -B.transpose() * info_matrix * e;
			// compute the contribution of this constraint to the linear system
			H_ii = A.transpose() * info_matrix * A;
			H_ij = A.transpose() * info_matrix * B;
			H_ji = B.transpose() * info_matrix * A;
			H_jj = B.transpose() * info_matrix * B;

			// construct the sparse matrix H and dense matrix b
			std::vector <Eigen::Triplet<DataType>> triplets;
			for( int i = 0; i < H_ii.rows(); i ++ ){
				for( int j = 0; j < H_ii.cols(); j ++ ){
					H_triplets.emplace_back( id_i * 3 + i, id_i * 3 + j, H_ii( i, j ) );			
				}
			}
			
			for( int i = 0; i < H_ij.rows(); i ++ ){
				for( int j = 0; j < H_ij.cols(); j ++ ){
					H_triplets.emplace_back( id_i * 3 + i, id_j * 3 + j, H_ij( i, j ) );
				}
			}
			
			for( int i = 0; i < H_ji.rows(); i ++ ){
				for( int j = 0; j < H_ji.cols(); j ++ ){
					H_triplets.emplace_back( id_j * 3 + i, id_i * 3 + j, H_ji( i, j ) );
				}
			}

			for( int i = 0; i < H_jj.rows(); i ++ ){
				for( int j = 0; j < H_jj.rows(); j ++ ){
					H_triplets.push_back( id_j * 3 + i, id_j * 3 + j, H_jj( i, j ) );
				}
			}
			
			//H.setFromTriplets( triplets.begin(), triplets.end() );
	
			b( id_i * 3 ) += b_i( 0 );
			b( id_i * 3 + 1 ) += b_i( 1 );
			b( id_i * 3 + 2 ) += b_i( 2 );
		
			b( id_j * 3 ) += b_j( 0 );
			b( id_j * 3 + 1 ) += b_j( 1 );
			b( id_j * 3 + 2 ) += b_j( 1 );
		}
	
		// keep the first node fixed
		H_triplets.emplace_back( 0, 0, 1 );
		H_triplets.emplace_back( 1, 1, 1 );
		H_triplets.emplace_back( 2, 2, 1 );
		
		H.setFromTriplets( H_triplets.begin(), H_triplets.end() );
	}

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
		
		Matrix3x3 zt_ij = v2t( z_ij );
		Matrix3x3 vt_i = v2t( v_i );
		Matrix3x3 vt_j = v2t( v_j );

		Matrix3x3 f_ij = vt_i.inverse() * vt_j;
	
		DataType theta_i = v_i[2];
		Vector2 t_i( v_i(0), v_i(1) );
		Vector2 t_j( v_j(0), v_j(1) );
	
		Vector2 dt_ij = t_j - t_i;

		DataType si = ::sin( theta_i );
		DataType ci = ::cos( theta_i );

		A << -ci, -si, -si * dt_ij[0] + ci * dt_ij[1],
		      si, -ci, -ci * dt_ij[0] - si * dt_ij[1],
		       0,   0,  1;
	
		B << ci, si, 0,
		    -si, ci, 0,
		      0,  0, 1;

		Matrix3x3 zt_ij_inv = zt_ij.inverse();
		e = t2v( zt_ij_inv * f_ij );
		e( 0, 2 ) = 0;
		e( 1, 2 ) = 0;
		
		A = zt_ij_inv * A;
		B = zt_ij_inv * B;
	}

private:
	const Matrix3x3 v2t( const Vector3 &v )
	{
		Matrix3x3 A;
	
		A << ::cos( v[2] ), -::sin( v[2] ), v[0],
		     ::sin( v[2] ),  ::sin( v[2] ), v[1],
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
	Eigen::SparseMatrix<DataType> H;	
	//Eigen::SparseMatrix<DataType> b;
	//Eigen::SparseMatrix<DataType> delta_x;
	//Eigen::SparseMatrix<DataType> x;
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
