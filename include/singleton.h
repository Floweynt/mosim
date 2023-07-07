#pragma once

template <class T>
class singleton
{
public:
    singleton& operator=(const singleton&) = delete;
    singleton& operator=(singleton&&) = delete;

    static T& get_instance()
    {
        if (!instance)
            instance = new inst;
        return *instance;
    }

protected:
    singleton() = default;

private:
    struct inst : public T
    {
        inst() : T() {}
    };

    static inline T* instance = nullptr;
};

