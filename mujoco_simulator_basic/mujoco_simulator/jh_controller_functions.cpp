#include "jh_controller.h"

MatrixXd mj2eigen(mjtNum *num, int nrow, int ncol) {
	MatrixXd res;
	res.resize(nrow, ncol);
	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			res(i, j) = num[i*ncol + j];
		}
	}
	return res;
}

mjtNum* eigen2mj(MatrixXd mat, mjData *d) {	
	Index matsize = mat.size();
	mjtNum *num = mj_stackAlloc(d, (int)matsize);
	for (int i = 0; i < mat.rows(); i++) {
		for (int j = 0; j < mat.cols(); j++) {
			 num[i*mat.cols() + j]=mat(i, j);
		}
	}
	return num;
}

MatrixXd Skm(Vector3d x) {
	Matrix3d Skew_temp1;
	Skew_temp1.resize(3, 3);
	Skew_temp1.setZero();
	Skew_temp1(0, 1) = -x(2);
	Skew_temp1(0, 2) = x(1);
	Skew_temp1(1, 0) = x(2);
	Skew_temp1(1, 2) = -x(0);
	Skew_temp1(2, 0) = -x(1);
	Skew_temp1(2, 1) = x(0);
	return Skew_temp1;
}

double Cubic(double rT, double rT_0, double rT_f, double rx_0, double rx_dot_0, double rx_f, double rx_dot_f) {
	double rx_t;
	if (rT<rT_0)
	{
		rx_t = rx_0;
	}
	else if (rT >= rT_0 && rT<rT_f)
	{
		rx_t = rx_0 + rx_dot_0*(rT - rT_0)
			+ (3 * (rx_f - rx_0) / ((rT_f - rT_0) * (rT_f - rT_0)) - 2 * rx_dot_0 / ((rT_f - rT_0) * (rT_f - rT_0)) - rx_dot_f / ((rT_f - rT_0) * (rT_f - rT_0)))*(rT - rT_0)*(rT - rT_0)
			+ (-2 * (rx_f - rx_0) / ((rT_f - rT_0) * (rT_f - rT_0) * (rT_f - rT_0)) + (rx_dot_0 + rx_dot_f) / ((rT_f - rT_0) * (rT_f - rT_0) * (rT_f - rT_0)))*(rT - rT_0)*(rT - rT_0)*(rT - rT_0);
	}
	else
	{
		rx_t = rx_f;
	}
	return (rx_t);
}

Matrix3d Rotation_xyz(double rad_x, double rad_y, double rad_z) {
	Matrix3d Rx, Ry, Rz;
	Rx << 1, 0, 0, 0, cos(rad_x), -sin(rad_x), 0, sin(rad_x), cos(rad_x);
	Ry << cos(rad_y), 0, sin(rad_y), 0, 1, 0, -sin(rad_y), 0, cos(rad_y);
	Rz << cos(rad_z), -sin(rad_z), 0, sin(rad_z), cos(rad_z), 0, 0, 0, 1;
	return Rx*Ry*Rz;
}


/*
bool isRotationMatrix(Matrix3d &R)
{
	Matrix3d Rt;
	transpose(R, Rt);
	Matrix3d shouldBeIdentity = Rt * R;
	Matrix3d I = Matrix3d::eye(3, 3, shouldBeIdentity.type());
	return  norm(I, shouldBeIdentity) < 1e-6;
}
*/


// Calculates rotation matrix to euler angles
// The result is the same as MATLAB except the order
// of the euler angles ( x and z are swapped ).
Vector3d rotationMatrixToEulerAngles(Matrix3d &R)
{
	double sy = sqrt(R(0, 0) * R(0, 0) + R(1, 0) * R(1, 0));

	bool singular = sy < 1e-6; // If

	double x, y, z;
	if (!singular)
	{
		x = atan2(R(2, 1), R(2, 2));
		y = atan2(-R(2, 0), sy);
		z = atan2(R(1, 0), R(0, 0));
	}
	else
	{
		x = atan2(-R(1, 2), R(1, 1));
		y = atan2(-R(2, 0), sy);
		z = 0;
	}
	return Vector3d(x, y, z);
}



void jh_controller::part::name(mjModel *m, mjData *d, char *part_name) {
	part_id = mj_name2id(m, mjOBJ_BODY, part_name);
	Jac_point.resize(6, m->nv);
	Jac_point.setZero();
	Jac.resize(6, m->nv);
	Jac.setZero();
	Jac_COM.resize(6, m->nv);
	Jac_COM.setZero();
	Jac_Contact.resize(6, m->nv);
	Jac_Contact.setZero();

	Rotm.resize(3, 3);
	Rotm.setZero();

	Mass = m->body_mass[part_id];

	m_int = m;
	d_int = d;

	inertia.setZero(3, 3);

	di.setZero(3, 3);

}

void jh_controller::part::refresh() {

	

	mjtNum* jac_r_temp = mj_stackAlloc(d_int, 3 * m_int->nv);
	mjtNum* jac_p_temp = mj_stackAlloc(d_int, 3 * m_int->nv);
	mj_jacBody(m_int, d_int, jac_p_temp,jac_r_temp,part_id);
	Jac.block(0, 0, 3, m_int->nv) = mj2eigen(jac_p_temp, 3, m_int->nv);
	Jac.block(3, 0, 3, m_int->nv) = mj2eigen(jac_r_temp, 3, m_int->nv);
	mj_jacBodyCom(m_int, d_int, jac_p_temp, jac_r_temp, part_id);
	Jac_COM.block(0, 0, 3, m_int->nv) = mj2eigen(jac_p_temp, 3, m_int->nv);
	Jac_COM.block(3, 0, 3, m_int->nv) = mj2eigen(jac_r_temp, 3, m_int->nv);
	for (int i = 0; i < 3; i++) {
		xpos(i) = d_int->xpos[part_id * 3 + i];
		xipos(i) = d_int->xipos[part_id * 3 + i];
	}
	

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Rotm(i, j) = d_int->xmat[part_id * 9 + i * 3 + j];
		}
	}
	
	Quaterniond q;
	q.w() = m_int->body_iquat[part_id * 4];
	q.x() = m_int->body_iquat[part_id * 4 + 1];
	q.y() = m_int->body_iquat[part_id * 4 + 2];
	q.z() = m_int->body_iquat[part_id * 4 + 3];
	i_rot = q.normalized().toRotationMatrix();
	


	

	di(0, 0) = m_int->body_inertia[part_id * 3];
	di(1, 1) = m_int->body_inertia[part_id * 3+1];
	di(2, 2) = m_int->body_inertia[part_id * 3+2];

	inertia = i_rot* di * i_rot.transpose();

	

	

	Jac = Jac*Rotation2g;
	Jac_COM = Jac_COM*Rotation2g;







}

void jh_controller::part::SetContact(Vector3d contact_vector) {
	Vector3d con_pos;
	mjtNum* jac_r_temp = mj_stackAlloc(d_int, 3 * m_int->nv);
	mjtNum* jac_p_temp = mj_stackAlloc(d_int, 3 * m_int->nv);
	con_pos = xpos + Rotm*contact_vector;

	mjtNum con_pos_mj[3];
	con_pos_mj[0] = con_pos[0];
	con_pos_mj[1] = con_pos[1];
	con_pos_mj[2] = con_pos[2];

	mj_jac(m_int, d_int, jac_p_temp, jac_r_temp, con_pos_mj, part_id);

	Jac_Contact.block(0, 0, 3, m_int->nv) = mj2eigen(jac_p_temp, 3, m_int->nv);
	Jac_Contact.block(3, 0, 3, m_int->nv) = mj2eigen(jac_r_temp, 3, m_int->nv);

	Jac_Contact = Jac_Contact*Rotation2g;


}
void jh_controller::part::SetContact(double x, double y,double z) {
	Vector3d con_pos,contact_vector;
	contact_vector << x, y, z;
	mjtNum* jac_r_temp = mj_stackAlloc(d_int, 3 * m_int->nv);
	mjtNum* jac_p_temp = mj_stackAlloc(d_int, 3 * m_int->nv);
	con_pos = xpos + Rotm*contact_vector;

	mjtNum con_pos_mj[3];
	con_pos_mj[0] = con_pos[0];
	con_pos_mj[1] = con_pos[1];
	con_pos_mj[2] = con_pos[2];

	mj_jac(m_int, d_int, jac_p_temp, jac_r_temp, con_pos_mj, part_id);

	Jac_Contact.block(0, 0, 3, m_int->nv) = mj2eigen(jac_p_temp, 3, m_int->nv);
	Jac_Contact.block(3, 0, 3, m_int->nv) = mj2eigen(jac_r_temp, 3, m_int->nv);

	Jac_Contact = Jac_Contact*Rotation2g;

}
void jh_controller::part::SetPoint(double x, double y, double z) {
	Vector3d con_pos, contact_vector;
	contact_vector << x, y, z;
	mjtNum* jac_r_temp = mj_stackAlloc(d_int, 3 * m_int->nv);
	mjtNum* jac_p_temp = mj_stackAlloc(d_int, 3 * m_int->nv);
	con_pos = xpos + Rotm*contact_vector;

	mjtNum con_pos_mj[3];
	con_pos_mj[0] = con_pos[0];
	con_pos_mj[1] = con_pos[1];
	con_pos_mj[2] = con_pos[2];

	mj_jac(m_int, d_int, jac_p_temp, jac_r_temp, con_pos_mj, part_id);

	Jac_point.block(0, 0, 3, m_int->nv) = mj2eigen(jac_p_temp, 3, m_int->nv);
	Jac_point.block(3, 0, 3, m_int->nv) = mj2eigen(jac_r_temp, 3, m_int->nv);

	Jac_point = Jac_point*Rotation2g;

}
void jh_controller::part::SetPoint(Vector3d point_vector) {

	Vector3d con_pos;
	mjtNum* jac_r_temp = mj_stackAlloc(d_int, 3 * m_int->nv);
	mjtNum* jac_p_temp = mj_stackAlloc(d_int, 3 * m_int->nv);
	con_pos = xpos + Rotm*point_vector;

	mjtNum con_pos_mj[3];
	con_pos_mj[0] = con_pos[0];
	con_pos_mj[1] = con_pos[1];
	con_pos_mj[2] = con_pos[2];

	mj_jac(m_int, d_int, jac_p_temp, jac_r_temp, con_pos_mj, part_id);

	Jac_point.block(0, 0, 3, m_int->nv) = mj2eigen(jac_p_temp, 3, m_int->nv);
	Jac_point.block(3, 0, 3, m_int->nv) = mj2eigen(jac_r_temp, 3, m_int->nv);

	Jac_point = Jac_point*Rotation2g;
}

jh_controller::part::part() {
}