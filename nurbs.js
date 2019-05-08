import {vec3, vec4} from './gl-matrix/index.js';

export default class NURBS {
    static makeHomogenousCoordinates(x, y, z, weight) {
        const homogenous = vec4.fromValues(x, y, z, 1.0);
        vec4.scale(homogenous, homogenous, weight);
        return homogenous;
    }

    static projectHomogenousCoordinates(homogenous) {
        const projected = vec3.create();
        vec3.scale(projected, homogenous, (homogenous[3] == 0.0) ? 1.0 : 1.0/homogenous[3]);
        return projected;
    }

    constructor() {
        this.degree = 0;
        this.controlPoints = [];
        this.knotSpans = [];
        this.cumulativeKnotSpans = [];
        this.multiplicity = [];
    }

    updateKnotSpans() {
        let parameter = 0, multiplicity = 0;
        this.cumulativeKnotSpans = [];
        this.multiplicity = [];
        for(let knotIndex = 0; knotIndex < this.knotSpans.length; ++knotIndex) {
            parameter += this.knotSpans[knotIndex];
            this.cumulativeKnotSpans.push(parameter);
            if(this.knotSpans[knotIndex] > 0) {
                for(let i = knotIndex-multiplicity; i < knotIndex; ++i)
                    this.multiplicity[i] = multiplicity;
                multiplicity = 0;
            }
            this.multiplicity.push(++multiplicity);
        }
        let knotIndex = this.knotSpans.length-1;
        if(this.knotSpans[knotIndex] == 0)
            for(let i = this.knotSpans.length-multiplicity; i < knotIndex; ++i)
                this.multiplicity[i] = multiplicity;
        this.maxKnotIndex = this.cumulativeKnotSpans.length-this.degree-1;
    }

    knotIndexAtParameter(parameter) {
        let low = 0, high = this.cumulativeKnotSpans.length;
        while(low < high) {
            const mid = (low+high)>>1;
            if(this.cumulativeKnotSpans[mid] <= parameter)
                low = mid+1;
            else
                high = mid;
        }
        return Math.max(this.degree, Math.min(low-1, this.maxKnotIndex-1));
    }

    basisWeightFunction(parameter, knotIndex, degree) {
        const knotBegin = this.cumulativeKnotSpans[knotIndex],
              knotEnd = this.cumulativeKnotSpans[knotIndex+degree],
              divisor = knotEnd-knotBegin;
        return (divisor == 0.0) ? 1.0 : (parameter-knotBegin)/divisor;
    }

    basisFunction(parameter, knotIndex, degree=this.degree) {
        if(degree == 0)
            return (parameter >= this.cumulativeKnotSpans[knotIndex] && parameter < this.cumulativeKnotSpans[knotIndex+1]) ? 1 : 0;
        return this.basisWeightFunction(parameter, knotIndex, degree)*this.basisFunction(parameter, knotIndex, degree-1)+
            (1-this.basisWeightFunction(parameter, knotIndex+1, degree))*this.basisFunction(parameter, knotIndex+1, degree-1);
    }

    homogenousCoordinatesAt(parameter) {
        const result = vec4.create();
        for(let knotIndex = 0; knotIndex < this.maxKnotIndex; ++knotIndex)
            vec4.scaleAndAdd(result, result, this.controlPoints[knotIndex], this.basisFunction(parameter, knotIndex));
        return result;
    }

    deBoor(parameter, knotIndex) {
        knotIndex = Math.max(this.degree, Math.min(knotIndex, this.controlPoints.length-1));
        const points = this.controlPoints.slice(knotIndex-this.degree, knotIndex+1).map(controlPoint => vec4.clone(controlPoint));
        for(let r = 0; r < this.degree; ++r)
            for(let j = 0; j < this.degree-r; ++j)
                vec4.lerp(points[j], points[j], points[j+1], this.basisWeightFunction(parameter, knotIndex-this.degree+j+r+1, this.degree-r));
        return points[0];
    }

    deBoorIntermediateSteps(parameter, knotIndex) {
        const results = [this.controlPoints.slice(knotIndex-this.degree, knotIndex+1).map(controlPoint => vec4.clone(controlPoint))];
        for(let r = 0; r < this.degree; ++r) {
            results.push([]);
            for(let j = 0; j < this.degree-r; ++j) {
                results[r+1][j] = vec4.create();
                vec4.lerp(results[r+1][j], results[r][j], results[r][j+1], this.basisWeightFunction(parameter, knotIndex-this.degree+j+r+1, this.degree-r));
            }
        }
        return results;
    }

    // Wolfgang Böhm Algorithm
    insertKnotAt(parameter) {
        let knotIndex = this.knotIndexAtParameter(parameter);
        const knotSpan = parameter-this.cumulativeKnotSpans[knotIndex],
              points = this.controlPoints.slice(knotIndex-this.degree, knotIndex+1).map(controlPoint => vec4.clone(controlPoint));
        for(let j = 0; j < this.degree; ++j)
            vec4.lerp(points[j], points[j], points[j+1], this.basisWeightFunction(parameter, knotIndex-this.degree+j+1, this.degree));
        points.length = this.degree;
        ++knotIndex;
        this.controlPoints.splice(knotIndex-this.degree, this.degree-1, ...points);
        this.knotSpans[knotIndex] -= knotSpan;
        this.knotSpans.splice(knotIndex, 0, knotSpan);
        this.updateKnotSpans();
    }
};
