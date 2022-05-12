#ifndef __GRAPH_OPTIMIZE_H
#define __GRAPH_OPTIMIZE_H

#include <vector>
#include <Eigen/Dense>
#include <cmath>


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
		

	GraphOptimize()
	{

	}

	~GraphOptimize()
	{

	}

	void getHessianDerived( const std::vector<Vector3> &v_poses,
				const std::vector<int> &from_ids,
				const std::vector<int> &to_ids,
				const std::vector<Vector3> &e_means )
	{

	}

	void linearFactors( const std::vector<Vector3> &v_poses, 
			    const std::vector<int> &from_ids,
			    const std::vector<int> &to_ids,
			    const std::vector<Vector3> &e_means,
			    const int index,
			    Matrix3x3 &A,
			    Matrix3x3 &B,	
			    Vector3 &b)
	{
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
		Vector2 t_i = v_i.head<2>();
		Vector2 t_j = v_j.head<2>();
	
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
		     ::sin( v[2] ),  ::sin( v[2] ), V[1],
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
};

}

#endif
