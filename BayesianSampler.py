import ultranest
import ultranest.stepsampler

def UltranestSampler(parameters,likelihood,prior,nsteps,live_points,max_ncalls):
    
    sampler = ultranest.ReactiveNestedSampler(parameters, likelihood, prior,log_dir='output',resume=True)
    sampler.stepsampler = ultranest.stepsampler.SliceSampler(
        nsteps=nsteps,
        generate_direction=ultranest.stepsampler.generate_mixture_random_direction,
    )

    result = sampler.run(live_points,max_ncalls)
    flat_samples = sampler.results['samples']
    return flat_samples