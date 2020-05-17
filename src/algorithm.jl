function optimize(F_ul::Function, # upper level objective function
                  f_ll::Function, # lower level objective function
                  bounds_ul::Array,
                  bounds_ll::Array,
                  method::Algorithm
                  )

      problem = Problem(F_ul,f_ll,bounds_ul,bounds_ll)
      engine = method.engine

      method.options.debug && @info("Initializing population...")
      method.status.initial_time = time()
      engine.initialize!(problem, engine, method.parameters, method.status, method.information, method.options)

      #####################################
      # common methods
      #####################################
      status = method.status
      information = method.information
      options     = method.options
      update_state! = engine.update_state!
      final_stage!  = engine.final_stage!
      ###################################

      ###################################
      # store convergence
      ###################################
      convergence = State[]
      if options.store_convergence
            st = deepcopy(status)
            push!(convergence, st)
      end
      
      method.options.debug && @info("Starting main loop...")
      options.debug && display(status)


      status.iteration = 0
      status.stop = status.stop || engine.stop_criteria(status, information, options)
      while !status.stop
            status.iteration += 1

            update_state!(problem,engine,method.parameters,method.status,method.information,method.options,status.iteration)
            
            options.debug && display(status)


            if options.store_convergence
                  st = deepcopy(status)
                  push!(convergence, st)
            end
            
            status.stop = status.stop || engine.stop_criteria(status, information, options)
      
      end

      status.convergence = convergence

      final_stage!(status, information, options)
      method.status.final_time = time()

      return status

end
