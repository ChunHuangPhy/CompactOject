import ultranest
import ultranest.stepsampler

def UltranestSampler(parameters,likelihood,prior,step,live_points,max_calls):
    """UltraNest based nested sampler by given likelihood prior, and parameters.
    
    Args:
        parameters (array): parameters array that want to be constrained.
        likelihood (array): theta as input. likelihood function defined by user.
        prior (array): cube as input, prior function defined by user. please check our test_inference.ipynb
        to check how to define likelihood and prior.
        step (int): as a step sampler, define this inference want to devided to how many steps.
        live_points (int): define how many live points will be used to explore the whole
        parameters space.
        max_ncalls (int): define after how many steps the sampler will stop work.
        
    Returns:
        flat_samples (array): equal weighted samples of whole posteior parameter space,
        this run will generate a dirctory as 'output', please check the run# folder, and the
        chain dirctory, there is a 'equal_weighted_samples' file, that is same with flat_samples here.
        It will be easier to check if you are using clusters to do this inference.
        
    """
    sampler = ultranest.ReactiveNestedSampler(parameters, likelihood, prior,log_dir='output')
    sampler.stepsampler = ultranest.stepsampler.SliceSampler(
        nsteps=step,
        generate_direction=ultranest.stepsampler.generate_mixture_random_direction,
        # adaptive_nsteps=False,
        # max_nsteps=400
    )

    result = sampler.run(min_num_live_points=live_points,max_ncalls= max_calls)
    flat_samples = sampler.results['samples']
    return flat_samples

def UltranestSamplerResume(parameters,likelihood,prior,nsteps,live_points,max_calls):
    """UltraNest based nested sampler by given likelihood prior, and parameters. (resume true verion
    could restart you run from your previous stopped results)
    
    Args:
        parameters (array): parameters array that want to be constrained.
        likelihood (array): theta as input. likelihood function defined by user.
        prior (array): cube as input, prior function defined by user. please check our test_inference.ipynb
        to check how to define likelihood and prior.
        step (int): as a step sampler, define this inference want to devided to how many steps.
        live_points (int): define how many live points will be used to explore the whole
        parameters space.
        max_ncalls (int): define after how many steps the sampler will stop work.
        
    Returns:
        flat_samples (array): equal weighted samples of whole posteior parameter space,
        this run will generate a dirctory as 'output', please check the run# folder, and the
        chain dirctory, there is a 'equal_weighted_samples' file, that is same with flat_samples here.
        It will be easier to check if you are using clusters to do this inference.
        
    """
    sampler = ultranest.ReactiveNestedSampler(parameters, likelihood, prior,log_dir='output',resume=True)
    sampler.stepsampler = ultranest.stepsampler.SliceSampler(
        nsteps=nsteps,
        generate_direction=ultranest.stepsampler.generate_mixture_random_direction,
        # adaptive_nsteps=False,
        # max_nsteps=400
    )

    result = sampler.run(min_num_live_points=live_points,max_ncalls=max_calls)
    flat_samples = sampler.results['samples']
    return flat_samples