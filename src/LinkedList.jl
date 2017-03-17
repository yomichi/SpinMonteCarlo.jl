type LinkedListNode{T}
    prev :: Union{LinkedListNode{T}, Void}
    next :: Union{LinkedListNode{T}, Void}
    val :: T
    function LinkedListNode{T}(val::T)
        new(nothing, nothing, val)
    end
end

type LinkedList{T}
    head :: Union{LinkedListNode{T}, Void}
    last :: Union{LinkedListNode{T}, Void}
end
LinkedList{T}(::Type{T}) = LinkedList{T}(nothing, nothing)

import Base: push!, pop!, shift!, unshift!, insert!, deleteat!
import Base: start, next, done

function push!{T}(ll::LinkedList{T}, val::T)
    lln = LinkedListNode{T}(val)
    if ll.last == nothing
        ll.head = ll.last = lln
    else
        lln.prev = ll.last
        ll.last.next = lln
        ll.last = lln
    end
    return lln
end

function pop!(ll::LinkedList)
    lln = ll.last
    if lln == ll.head
        ll.head = ll.last = nothing
    else
        lln.prev.next = nothing
        ll.last = lln.prev
    end
    return lln
end

function unshift!{T}(ll::LinkedList{T}, val::T)
    lln = LinkedListNode{T}(val)
    if ll.head == nothing
        ll.head = ll.last = lln
    else
        lln.next = ll.head
        ll.head.prev = lln
        ll.head = lln
    end
    return lln
end

function shift!(ll::LinkedList)
    lln = ll.head
    if lln == ll.last
        ll.head = ll.last = nothing
    else
        lln.next.prev = nothing
        ll.head = lln.next
    end
    return lln
end

function insert!{T}(ll::LinkedList{T}, next::LinkedListNode{T}, val::T)
    if next == nothing
        return push!(ll,val)
    end
    if next == ll.head
        return unshift!(ll, val)
    end
    lln = LinkedListNode(val)
    P = next.prev
    lln.prev = P
    lln.next = next
    P.next = lln
    next.prev = lln
end

function deleteat!(ll::LinkedList, lln::LinkedListNode)
    lln.prev.next = lln.next
    lln.next.prev = lln.prev
end

start(ll::LinkedList) = ll.head
next(ll::LinkedList, lln::LinkedListNode) =  lln.val, lln.next
done(ll::LinkedList, lln) = lln == nothing
