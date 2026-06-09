# Why model-selection risk is the load-bearing idea

A simulation is only as trustworthy as its weakest, least-validated input. That sentence is the whole project. Everything else in Hypnos is plumbing in service of making it true by default.

Here is the problem the plumbing solves.

## Same drug, same dose, two answers that disagree by 136%

Give a 72-year-old, 60 kg woman a 2 mg/kg bolus of propofol followed by a 6 mg/kg/h infusion. Ask two of the most widely deployed pharmacokinetic models what her effect-site concentration will be, the concentration at the brain that actually drives loss of consciousness.

The Marsh model says the peak is about 3.7 µg/mL. The Schnider model says about 8.2 µg/mL. Same patient. Same syringe. The two answers differ by 5.6 µg/mL, which is 136% of their average.

That is not a rounding error. In propofol terms it is the difference between comfortably asleep and dangerously deep. The gap comes mostly from one parameter, the plasma-to-effect-site equilibration rate constant written *k*e0, which the two models simply estimated differently. Marsh uses 0.26 per minute, Schnider 0.456. Neither is wrong. They were fit to different data with different methods, and the field has never agreed on one.

Now here is the uncomfortable part. If you run a single-model target-controlled infusion simulator, you see one curve. Confident, smooth, plotted to three significant figures. The disagreement is invisible. The simulator does not lie to you, exactly. It just declines to mention that the model down the hall would have told you something 136% different.

The field has said this out loud, in print: the availability of multiple PK/PD models for a single drug increases the risk of invalid model selection by the user. Hypnos takes that sentence literally and builds the missing instrument around it.

## Make the disagreement a measurement, not a footnote

The headline feature of Hypnos is a comparison, not a prediction. Pick a virtual patient and a dose. Hypnos overlays the predicted curves from every eligible model, greys out the ones the patient falls outside of, and reports the divergence as a number.

This reframes the question. The bad question is "what should I give this patient," which is a dosing decision and a regulated medical act. The good question is "how much do published models disagree here," which is a research question with a quantitative answer. Hypnos only ever answers the second one. That boundary is not a disclaimer bolted on at the end. It is enforced by construction: there is no inverse control anywhere in the code, no function that computes the infusion required to reach a target. The arithmetic that turns a dataset into a dosing device is simply absent.

## The envelope is where models go to break

Every population model was fit to a particular set of patients: a range of ages, weights, heights. Outside that range the model is not predicting, it is extrapolating, and extrapolation is where the quiet failures live.

Hypnos makes the derivation range a first-class, machine-readable field called the applicability envelope, and the simulator enforces it. Step outside, and the result is automatically tiered down to D, the lowest confidence tier, with a warning attached. You cannot accidentally get an A-looking number from an out-of-envelope guess. The honesty is structural, not editorial.

The sharpest example is the Schnider model in an obese patient. Schnider scales clearance using lean body mass computed by the James equation. The James equation has a mathematical wart: as weight climbs at a fixed height, computed lean body mass rises, peaks, and then *falls*. Combine that downturn with the negative coefficient on lean body mass in Schnider's clearance equation, and the model predicts that a heavier patient clears propofol faster, which is not physiology, it is an artifact. For a 140 kg patient the curve spikes to a non-physical place. Hypnos knows this, because the failure mode is encoded as a cited, machine-evaluable predicate (`bmi > 42`), and the obese patient's Schnider curve is greyed out and labeled rather than quietly plotted.

Children are the other classic trap. Most propofol models were built on adults. Run one on a six-year-old and you are not doing pediatric anesthesia, you are doing wishful arithmetic. Hypnos names that specifically: a sub-range patient in an adult model is flagged as a pediatric extrapolation and tiered to D, while the one model actually built for children, Paedfusor, is the only curve left in color. Same logic runs the other direction, and for the elderly. The label is not "out of envelope." The label is "an adult model used in a child is not predictive."

## Worst input wins

A real depth-of-anesthesia simulation is a stack: a PK model feeds a *k*e0 link, which feeds a pharmacodynamic effect model, sometimes plus a drug-drug interaction surface. Each layer has its own confidence tier.

Hypnos propagates the worst one. A Tier-A pharmacokinetic model driven through a Tier-C effect model yields a Tier-C answer, and the dashboard, the API, and the provenance baked into every export all say C. You do not get to average your confidence. A chain is as strong as its weakest link, and a simulation is as trustworthy as its least-validated component. Making that arithmetic automatic is most of the point.

## Honest about itself, out loud

Here is the line that earns the project its self-respect: humans verify, and the language model does not promote.

A model record moves from `unverified` to `verified` only when a person opens the source PDF and confirms, field by field, every structural parameter and, more importantly, every covariate equation. The covariate equations are where transcription errors hide, because a single flipped sign or a wrong lean-body-mass formula propagates silently into every downstream number. No automated agent gets to mark something verified on its own authority.

So the dataset reports its own status without flinching. Right now all fifteen models read `unverified`, and `hypnos status` will tell you exactly that: zero of fifteen verified, here are the highest-leverage ones to check first, here is the field-by-field checklist for each. Even the models that do run get no green checkmark. The general-purpose Eleveld propofol kernel was transcribed from the published equations, cross-checked against an open implementation, and validated to reproduce its reference patient to the decimal, and it still reads `unverified`, because reproducing one reference point is not the same as a human confirming every covariate equation against the source PDF. That is the whole discipline in one example.

A few models go further and refuse to run at all. The rocuronium pharmacokinetic model is curated and carries a verified citation, but its compartmental parameters cannot be reconciled cleanly across sources, so `simulate()` declines to evaluate it rather than risk shipping numbers the original paper would not recognize. Refusing to compute is, occasionally, the most honest thing a tool can do.

## Where pharmacometrics gives the tiers teeth

Confidence tiers could have been pure editorial opinion. They are not, entirely, because anesthesia comes with standard, quantitative accuracy metrics. Two matter most. Median performance error (MDPE) measures bias, whether a model systematically runs high or low. Median absolute performance error (MDAPE) measures inaccuracy, how far off it is regardless of direction. The Eleveld model, for instance, reports a roughly minus 27% bias in older adults, meaning it tends to underpredict their plasma concentrations, while keeping absolute inaccuracy under 30% across groups. Numbers like these let tier assignment be partly arithmetic instead of purely a matter of taste, and Hypnos records them, with citations, in each model's `predictive_performance`.

## The case against, stated plainly

Hypnos is infrastructure, and infrastructure has costs and limits worth naming.

It is not exhaustive and is not trying to be. The scope has a declared ceiling: published population models for anesthetic drugs in humans, full stop. Growth means adding enumerated models inside that envelope, not sprawling outward. Disease-state pharmacokinetics, novel agents without a peer-reviewed model, veterinary work, proprietary device algorithms: all explicitly out, by design, not by oversight.

It is not a simulator's replacement and not an automated researcher. It does forward simulation for comparison and export validation, nothing fancier.

And the honesty stance has a real price. Refusing to run the rocuronium model means a user who wants neuromuscular-block simulation today cannot get it from Hypnos today. Greying out Schnider for an obese patient means someone who wanted that number anyway has to go elsewhere to be misled. Marking everything `unverified`, even the kernels that reproduce their reference patient to the decimal, is less impressive than a green checkmark would be. These are deliberate tradeoffs. The project would rather under-claim and be right than over-claim and be trusted incorrectly, because the drugs at its core are lethal at the wrong dose and the failure mode of a confident wrong number is not embarrassment, it is harm.

## The opportunity

What this unlocks is a shared, honest substrate. The field already has open raw data on one side and individual published models on the other, with nothing curated, tiered, and machine-readable sitting between them. Hypnos is that middle layer: the thing a TCI researcher, a closed-loop control study, or a machine-learning-on-physiological-data project can build on without each re-typing the Schnider table from a 1998 PDF and re-introducing the same transcription bug the last three papers had.

Make model-selection uncertainty a first-class, machine-readable field, and it stops being a footnote that careful readers notice and careless ones miss. It becomes a number you can see, sort, enforce, and export.

That is the load-bearing idea. A simulation is only as trustworthy as its weakest, least-validated input, so make that fact the most visible thing in the room.
