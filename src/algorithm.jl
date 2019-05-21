function optimize(F_ul::Function, # upper level objective function
                  f_ll::Function, # lower level objective function
                  bounds_ul::Array,
                  bounds_ll::Array,
                  method::Algorithm,
                  information::Information = Information()
                  )

      status = method.initial_state
      information = method.information
      options = method.options
      final_stage! = method.final_stage!
      
      t::Int = 0
      while !method.stop_criteria(status, information, options)
            t += 1

            update_state!(F_ul,f_ll,bounds_ul,bounds_ll,status,information,method.options,t)
            
            debug && display(status)
            
      end

      final_stage!(status, information, method.options)

      return status

end
