function optimize(F_ul::Function, # upper level objective function
                  f_ll::Function, # lower level objective function
                  bounds_ul::Array,
                  bounds_ll::Array,
                  method::Algorithm
                  )

      problem = Problem(F_ul,f_ll,bounds_ul,bounds_ll)
      engine = method.engine
      engine.initialize!(problem, method.parameters, method.status, method.information, method.options)

      #####################################
      # common methods
      #####################################
      status = method.status
      information = method.information
      options     = method.options
      update_state! = engine.update_state!
      final_stage!  = engine.final_stage!
      ###################################
      
      status.iteration = 0
      while !engine.stop_criteria(status, information, options)
            status.iteration += 1

            update_state!(problem,engine,method.parameters,method.status,method.information,method.options,status.iteration)
            
            options.debug && display(status)


            if options.store_convergence
                  st = deepcopy(status)
                  empty!(st.convergence)
                  push!(status.convergence, st)
            end
            
      end

      final_stage!(status, information, options)

      return status

end
