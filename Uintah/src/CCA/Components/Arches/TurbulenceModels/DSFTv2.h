#ifndef Uintah_Component_Arches_DSFTv2_h
#define Uintah_Component_Arches_DSFTv2_h

#include <CCA/Components/Arches/Task/TaskInterface.h>
#include <CCA/Components/Arches/TurbulenceModels/DynamicSmagorinskyHelper.h>

namespace Uintah{

  class DSFTv2 : public TaskInterface {

public:

    DSFTv2( std::string task_name, int matl_index );
    ~DSFTv2();

    void problemSetup( ProblemSpecP& db );

    void register_initialize( std::vector<ArchesFieldContainer::VariableInformation>& variable_registry , const bool packed_tasks);

    void register_timestep_init( std::vector<ArchesFieldContainer::VariableInformation>& variable_registry , const bool packed_tasks);

    void register_timestep_eval( std::vector<ArchesFieldContainer::VariableInformation>& variable_registry, const int time_substep , const bool packed_tasks);

    void register_compute_bcs( std::vector<ArchesFieldContainer::VariableInformation>& variable_registry, const int time_substep , const bool packed_tasks){}

    void compute_bcs( const Patch* patch, ArchesTaskInfoManager* tsk_info ){}

    void initialize( const Patch* patch, ArchesTaskInfoManager* tsk_info );

    void timestep_init( const Patch* patch, ArchesTaskInfoManager* tsk_info );

    void eval( const Patch* patch, ArchesTaskInfoManager* tsk_info );

    void create_local_labels();

    //Build instructions for this (DSFTv2) class.
    class Builder : public TaskInterface::TaskBuilder {

      public:

      Builder( std::string task_name, int matl_index ) : m_task_name(task_name), m_matl_index(matl_index){}
      ~Builder(){}

      DSFTv2* build()
      { return scinew DSFTv2( m_task_name, m_matl_index ); }

      private:

      std::string m_task_name;
      int m_matl_index;
    };

private:
    std::string m_u_vel_name;
    std::string m_v_vel_name;
    std::string m_w_vel_name;
    std::string m_density_name;
    std::string m_volFraction_name;

    std::string m_cc_u_vel_name;
    std::string m_cc_v_vel_name;
    std::string m_cc_w_vel_name;

    std::string m_rhou_vel_name;
    std::string m_rhov_vel_name;
    std::string m_rhow_vel_name;
    std::string m_IsI_name;
    std::string m_ref_density_name;
    std::string m_cell_type_name;
    //int Type_filter ;
    bool m_create_labels_IsI_t_viscosity{true};
    Uintah::ArchesCore::FILTER Type_filter;
    Uintah::ArchesCore::TestFilter m_Filter;
  };
}
#endif
