
//  Variables for WH :
//  InputVarNames = jet0_times_jet1 jet0_times_lep jet1_times_lep jet0_times_nu jet1_times_nu lep_times_nu gammaWH_zaxis angle_Bb_beam angle_plane_Bb_beam_plane_Bb_lep
//  
//  Variables for ZH
//  InputVarNames = jet0_times_jet1 jet0_times_lep jet1_times_lep jet0_times_lep2 jet1_times_lep2 lep_times_lep2 gammaZH_zaxis angle_Bb_beam angle_plane_Bb_beam_plane_Bb_lep
//  
//  differences wrt. WH :
//  jet0_times_nu ->  jet0_times_lep2
//  jet1_times_nu ->  jet1_times_lep2
//  lep_times_nu  ->  lep_times_lep2
//  gammaWH_zaxis ->  gammaZH_zaxis


#include <TLorentzVector.h>
#include <TVector2.h>

class LorentzVars {
private:

  TLorentzVector m_jet0;
  TLorentzVector m_jet1;
  TLorentzVector m_lep;
  TLorentzVector m_lep2;
  TVector2 m_met;

public:

  LorentzVars() {
  }

  void setVectorsWlvH(TLorentzVector jet0, TLorentzVector jet1, TLorentzVector lep, TVector2 met) {
    // jet0   highest pt (b-tagged) jet
    // jet1   second highest pt (b-tagged) jet
    // lep    the only candidate electron or muon
    // met    missing transverse energy
    m_jet0 = jet0;
    m_jet1 = jet1;
    m_lep = lep;
    m_met = met;
  }

  void setVectorsZllH(TLorentzVector jet0, TLorentzVector jet1, TLorentzVector lep, TLorentzVector lep2) {
    // jet0   highest pt (b-tagged) jet
    // jet1   second highest pt (b-tagged) jet
    // lep    the highest pt candidate electron or muon
    // lep2   the second highest pt candidate electron or muon
    m_jet0 = jet0;
    m_jet1 = jet1;
    m_lep = lep;
    m_lep2 = lep2;
  }

  // ============ WH variables ============

  float jet0_times_jet1() const {
    // 4-vector scalar product
    double val = m_jet0 * m_jet1;
    return val / (7.E9 + val);
  }

  float jet0_times_lep() const {
    // 4-vector scalar product
    double val = m_jet0 * m_lep;
    return val / (1.13E10 + val);
  }

  float jet1_times_lep() const {
    // 4-vector scalar product
    double val = m_jet1 * m_lep;
    return val / (5.2E9 + val);
  }

  TLorentzVector estimated_nu() const {
    // rough estimation of neutrino four vector

    // first set only transversal component from MET
    TLorentzVector nu(m_met.Px(), m_met.Py(), 0, 0);

    // then add z component by assuming nu takes on average 25% of momentum in decay
    TLorentzVector sum_Bb_lep(m_jet0 + m_jet1 + m_lep);
    nu += TLorentzVector(0, 0, 0.25 / 0.75 * sum_Bb_lep.Pz(), 0);

    // now add new energy calculated from full vector
    nu += TLorentzVector(0, 0, 0, nu.P());
    return nu;
  }

  float jet0_times_nu() const {
    // 4-vector scalar product
    double val = m_jet0 * estimated_nu();
    return val / (9.4E9 + val);
  }

  float jet1_times_nu() const {
    // 4-vector scalar product
    double val = m_jet1 * estimated_nu();
    return val / (4.5E9 + val);
  }

  float lep_times_nu() const {
    // 4-vector scalar product
    double val = m_lep * estimated_nu();
    return val / (4.6E9 + val);
  }

  float gammaWH_zaxis() const {
    // z-boost
    TLorentzVector sum_WH(m_jet0 + m_jet1 + m_lep + estimated_nu());
    return fabs(sum_WH.Pz()) / sum_WH.M();
  }

  float angle_Bb_beam() const {
    // angle
    TVector3 beam(0., 0., 1.);
    TVector3 vecBb((m_jet0 + m_jet1).Vect());
    return vecBb.Angle(beam);
  }

  float angle_plane_Bb_beam_plane_Bb_lep() const {
    // angle

    // get 3-vector using the 3-spatial components
    TVector3 beam(0., 0., 1.);
    TVector3 vecBb((m_jet0 + m_jet1).Vect());
    TVector3 Bb_cross_beam(vecBb.Cross(beam));

    // then 3-vector cross product with only or highest pt electron or muon candidate
    TVector3 Bb_cross_lep(vecBb.Cross(m_lep.Vect()));
    double angle = Bb_cross_beam.Angle(Bb_cross_lep);
    return (Bb_cross_beam.Cross(Bb_cross_lep)).Dot(vecBb) > 0 ? angle : -angle;
  }


  // ============ ZH variables ============
  // the following WH variables are replaced :
  // jet0_times_nu ->  jet0_times_lep2
  // jet1_times_nu ->  jet1_times_lep2
  // lep_times_nu  ->  lep_times_lep2
  // gammaWH_zaxis ->  gammaZH_zaxis

  float jet0_times_lep2() const {
    // 4-vector scalar product
    double val = m_jet0 * m_lep2;
    return val / (9.4E9 + val);
  }

  float jet1_times_lep2() const {
    // 4-vector scalar product
    double val = m_jet1 * m_lep2;
    return val / (4.5E9 + val);
  }

  float lep_times_lep2() const {
    // 4-vector scalar product
    // m_wCandidate.lepton()  the highest pt candidate electron or muon
    // m_wCandidate.lepton2() the second highest pt candidate electron or muon
    double val = m_lep * m_lep2;
    return val / (4.6E9 + val);
  }

  float gammaZH_zaxis() const {
    // z-boost
    TLorentzVector sum_ZH(m_jet0 + m_jet1 + m_lep + m_lep2);
    return fabs(sum_ZH.Pz()) / sum_ZH.M();
  }

};
