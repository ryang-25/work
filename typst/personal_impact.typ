#set page("us-letter")

#text(rgb("#CD9575"))[= Improving Computer Vision Accuracy]

#line(length: 100%)

25 July 2024

Roland Yang


Although it may seem like a recent innovation, computer vision has a rich
history spanning back to the 1960s. Automating the process of recognizing
objects with human-level accuracy has enabled large swaths of tasks to be done
autonomously and with minimal oversight. Now, it seems as though computer
vision is prolific in our lives. Our cars are able to seemingly "see" objects
and drive independently, our smartphones are able to recognize the subject of a
photograph, and rare medical diseases can be detected from films with
remarkable accuracy.

Despite such widespread deployment of the technology, computer vision maintains
countless deficiencies. While these models recognize images by picking up on
small distinguishing details, their accuracy can plummet if the image data is
not perfectly formatted, and in the real world it often isn't. Photos may be
taken in bad lighting, at different angles, or in hazy weather conditions.
Computer vision models are expected to accommodate all of these, maintaining
their accuracy as a human would. In reality, computer models function contrary
to our expectations. They function poorly when inputs are presented in an
unexpected format, and they can be tricked rather trivially. Perform the right
set of manipulations to an image and a model might confidently predict that a
picture of a horse contains a dog. The danger of such vulnerabilities can be
great when image classification is used in safety-critical scenarios and is why
remediation is an ongoing area of research.

When we deploy computer vision in the real-world, we may want to deploy defenses
to render our models more resilient to image manipulations. Although there
exist a plethora of ``training recipes'' that do just that, there is a want to
be able to compare them, since not are equally effective and make different
tradeoffs. Some may sacrifice the accuracy on unmanipulated inputs to boost
accuracy on manipulated inputs, and what users desire is contingent on the
exact problem domain that the model is used in. Comparing on accuracy alone is
insufficient, since there isn't enough data to capture how each defense
performs. More metrics help greatly, both in evaluating the efficacy of attacks
and informing the design of better defenses.

Our contribution is then simple: we implement several metrics for a number of
defenses to evaluate how they compare. Our choice of which metrics to evaluate
is informed by research, and our chosen defenses are similarly well regarded as
effective. We run multiple trials on different architectures and datasets to demonstrate
that defenses generalize well to a variety of models. We find that our
results are an improvement on previously published work; furthermore, we see
that these techniques have wider applications than initially thought, which may
be of great interest for future work to build upon. Our work will be published
in #link("https://arxiv.org")[arXiv] soon.

On this topic, I recently had a look at the slides of "Guaranteed Safe AI and
Robustness" from the Alignment Workshop in Vienna where Nicholas Carlini, one
of the most foremost researchers in this field, mentioned that there have been
over 9,000 papers published on this topic within the last 10 years. We have yet
to realize a satisfactory solution to this problem of computer vision being
vulnerable in a scenario where attackers are allowed unfettered access to the
model. Although these attacks are now being shown to be practically feasible, I
believe that there should be a greater emphasis still on helping models to
classify accurately under a wider range of conditions, since that work is most
pertinent in the real world, because, let's face it: how often will your Tesla
be tricked by a street sign that someone wanted to make appear as a stop sign?
Likely never. While it is a worthwhile goal, I think that your car being able
to drive safely in snow or fog is of a larger concern. For those that desire a
bit of improvement in both categories, our work shows that existing techniques
have the potential to address both, albeit with reduced effectiveness in each.
By converging these two fields in our work, it sets the framework for more
opportunities for experimentation and further study.