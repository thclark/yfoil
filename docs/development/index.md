# Development

!!! note "Stub — to be written"

    How yFoil is built, kept faithful to XFOIL, and tested. The [conventions](conventions/naming.md)
    and the [testing](testing/README.md) reports are in the sidebar.

## Use of AI

As will be very clear to anyone familiar with the use of coding agents, yFoil's development is **heavily** dependent on the use of agentic AI tools. 

There's simply no other practical way of marshalling the sheer amount of data required for validation and unit test extraction... covering every line,
and (much harder!) every scenario of logical branch execution available in a program as complex as XFOIL, would be spectacularly hard. Without AI it
would take years, possibly decades, of software development effort to produce work of this scope and quality.

Agentic AI comes under a lot of fire for 'vibe coding' of poor quality, however, this is not that: Reconstruction like this is the perfect use 
case - agents thrive on constraints, and there's no better constraint than "it must produce this number from a translation of that original source code".

Our response to the phrase "vibe coding" is that, if this is vibe coding, then it's an extremely demanding, militant vibe... desiring only the highest quality!

### Where it lets us down

Document generation. If you're reading a coherent part of the documentation, that's because we've whipped outputs into the right
places and carefully worded it. If you're reading something that's a bit more like words splattered on a page - it's an AI placeholder. We'll slowly get
to it until the entire set of docs is well polished.

### Agents used

Anthropic's Claude Code harness was used exclusively, with:

- Opus < 5.4 model for early data structures and iterations (it didn't cope well)
- Fable 5.x model delivering the core implementation (very capable but costly, exhausted a 20x max plan several times over a couple of months)
- Opus 5.5 model reworking validation studies and tests into a coherent structure (used since its release, our experience is that it is at least as capable as Fable, but quicker and much cheaper) 

Token budget is difficult to estimate as we weren't recording as we went, but the majority of the project was achieved working ~40% human time over ~3 months, on a Claud Max 20x subscription (175GBP/month).
