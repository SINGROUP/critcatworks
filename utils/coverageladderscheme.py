# scheme for coverage ladder workflow
# pure layout of logic



        # ladder variables: l = number of minimum energies selected
        #                   d = maximum depth of branches (convergence criterion)
        #                   n0 = starting number of adsorbates
        #                   structure0 = starting structure with adsorbates
        


        # DAG dictionary of jobs, simulation_ids

        # dictionary of jobs per number of adsorbates (key = n, value = simulation_id)

        # tracing variable how many adsorbates on the surface

        # adsorption energy computation

        # if root calculation
        # decision to go up or down the ladder
        # based on adsorption energy
        # based on depth

        # -------------------------------#

        # similarity of occupied / unoccupied sites
        # array of new structures with one adsorbate added / removed

        # spawning of new jobs

        # --------------------------------#

        # select l minimum energy structures (l is fixed parameter) 
        # comparing among energies of same number of adsorbates n
        # if all in branch too high, continue on another branch
        # otherwise, select as new root

        # 




def main(l = 2, d = 4, max_steps = 100,):
    for step_number in range(max_steps):
        pass
        run_dft_jobs()

        direction, structure = decide_next_branch()

        if direction == "up":
            add_adsorbate()
        else:
            remove_adsorbate()

        




    return


if __name__ == '__main__':
    main()
